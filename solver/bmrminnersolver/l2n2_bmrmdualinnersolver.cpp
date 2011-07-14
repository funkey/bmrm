/* Copyright (c) 2009, NICTA
 * All rights reserved.
 *
 * The contents of this file are subject to the Mozilla Public License Version
 * 1.1 (the "License"); you may not use this file except in compliance with
 * the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS IS" basis,
 * WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License
 * for the specific language governing rights and limitations under the
 * License.
 *
 * Authors      : Choon Hui Teo (ChoonHui.Teo@anu.edu.au)
 * Created      : 20/11/2007
 * Last Updated :
 */

#ifndef _L2N2_BMRMINNERSOLVER_CPP_
#define _L2N2_BMRMINNERSOLVER_CPP_

#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>
#include "bmrmexception.hpp"
#include "l2n2_bmrmdualinnersolver.hpp"
#include "configuration.hpp"

// when density of gradient exceeds this value, store the gradient in dense format
#define SPARSE_TO_DENSE_THRESHOLD  0.6
#define AGG_GRAD_TIME_STAMP 99999

using namespace std;


// in general, BMRM needs more number of iterations to converge if we remove all IDLE gradients at once

CL2N2_BMRMDualInnerSolver::CL2N2_BMRMDualInnerSolver(double lambda)
   : CBMRMInnerSolver(lambda),
     prevDim(0),
     f(0),
     Q(0),
     a(0),
     b(0),
     l(0),
     u(0),
     tol(1e-6),
     removeAllIdleGrad(false),
     QPScale(1.0),
     aggGradIdx(-1),
     useDiagQuadForm(false),
     B(0),
     invB(0)
{
   // inherited from CInnerSolver
   iter = 0;
   dim = 0;
   numOfConstraint = 1;
   
   // set (default) private member values
   gradIdleAge = 10;
   maxGradSetSize = 2048;  // max num of gradients
   
   // synchronize members with configuration file
   Configuration &config = Configuration::GetInstance();
   
   if(config.IsSet("L2N2_BMRMDualInnerSolver.maxGradSetSize"))
      maxGradSetSize = config.GetInt("L2N2_BMRMDualInnerSolver.maxGradSetSize");
   
   if(config.IsSet("L2N2_BMRMDualInnerSolver.gradIdleAge")) 
   {
      gradIdleAge = config.GetInt("L2N2_BMRMDualInnerSolver.gradIdleAge");
      gradIdleAge = max(gradIdleAge, 2);
   }
   
   if(config.IsSet("L2N2_BMRMDualInnerSolver.removeAllIdleGradients"))
      removeAllIdleGrad = config.GetBool("L2N2_BMRMDualInnerSolver.removeAllIdleGradients");
   
   if(config.IsSet("L2N2_BMRMDualInnerSolver.diagonalPDMatrix"))
   {
      string diagPDMatrixFn = "";
      diagPDMatrixFn = config.GetString("L2N2_BMRMDualInnerSolver.diagonalPDMatrix");
      vector<double> tmpB;
      
      FILE *ifp = fopen(diagPDMatrixFn.c_str(),"r");
      if(ifp == NULL)
         throw CBMRMException("Cannot open diagnoalPDMatrix\n","CL2N2_BMRMDualInnerSolver::CL2N2_BMRMDualInnerSolver()");
      
      double tmp = 0;
      while(!feof(ifp))
      {
         fscanf(ifp,"%lf\n",&tmp);
         
         if(tmp <= 0)
            throw CBMRMException("Diagonal PD Matrix should be strictly positive\n","CL2N2_BMRMDualInnerSolver::CL2N2_BMRMDualInnerSolver()");
         tmpB.push_back(tmp);
      }
      fclose(ifp);
      useDiagQuadForm = true;
      
      B = new TheMatrix(1,tmpB.size(),SML::DENSE);
      invB = new TheMatrix(1,tmpB.size(),SML::DENSE);
      
      for(size_t i=0; i<tmpB.size(); i++)
      {
         B->Set(i,tmpB[i]);
         invB->Set(i, 1.0/tmpB[i]);
      }
      
      if(verbosity > 1)
      {
         double maxi = -SML::INFTY;
         double mini = SML::INFTY;
         for(size_t i=0; i<tmpB.size(); i++)
         { 
            double tmp = 0;
            B->Get(i,tmp);
            maxi = max(tmp,maxi);
            mini = min(tmp,mini);
         }
         cout << "max(B) = " << maxi << endl;
         cout << "min(B) = " << mini << endl;
      }
   }
   
   // pre-allocate memory for gradient, offset and active and enter time-stamps
   gradientSet.resize(maxGradSetSize,0);                  
   offsetSet.resize(maxGradSetSize,0);                  
   activeTimeStamp.resize(maxGradSetSize,0);  
   enterTimeStamp.resize(maxGradSetSize,0);   
   
   // allocate maximum memory needed
   x = (double*)calloc(maxGradSetSize, sizeof(double));
   Q = (double*)calloc(maxGradSetSize*maxGradSetSize, sizeof(double));
   f = (double*)calloc(maxGradSetSize, sizeof(double));
   l = (double*)calloc(maxGradSetSize, sizeof(double));
   u = (double*)calloc(maxGradSetSize, sizeof(double));
   a = (double*)calloc(maxGradSetSize, sizeof(double));
   b = new double(1.0); 
   
   // initializations
   Reset();
}


CL2N2_BMRMDualInnerSolver::~CL2N2_BMRMDualInnerSolver()
{
   // free for parent
   if(x) free(x);
   if(Q) free(Q);
   if(f) free(f);
   if(l) free(l);
   if(u) free(u);
   if(a) free(a);
   if(b) delete b;
   
   for(int i=0; i < dim; i++)
      if(gradientSet[i])
         delete gradientSet[i];   
   
   if(B) delete B;
   if(invB) delete invB;
}


void CL2N2_BMRMDualInnerSolver::Reset()
{
   // QP mem and gradient related mem
   for(int i=0; i<maxGradSetSize; i++) 
   {
      x[i] = 0;
      f[i] = 0;                
      l[i] = 0;
      u[i] = 1.0;
      a[i] = 1.0;      
      
      if(gradientSet[i]) 
      {
         delete gradientSet[i];
         gradientSet[i] = 0;
      }
      offsetSet[i] = 0.0;
      activeTimeStamp[i] = 0;
      enterTimeStamp[i] = 0;
   }                
   memset(Q, 0, sizeof(double)*maxGradSetSize*maxGradSetSize); 
   
   // variables
   iter = 0;
   dim = 0;
}


void CL2N2_BMRMDualInnerSolver::ChangeLambda(double newLambda)
{
   // rescale QP matrix and vector
   double change = bmrmLambda/newLambda;
   for(int i=0; i<dim; i++)
   {
      for(int j=0; j<dim; j++)
         Q[i*dim+j] *= change;
   }   
   bmrmLambda = newLambda;
}

void CL2N2_BMRMDualInnerSolver::Update(TheMatrix& a, double b)
{
   iter++;
   int idx = 0;
   prevDim = dim;
   
   if(removeAllIdleGrad)
      idx = RemoveAllIdleGradients();
   else
      idx = RemoveLaziestIdleGradient();
   
   
   if(dim > maxGradSetSize)
   {
      dim = maxGradSetSize;
      idx = AggregateGradients();
   }
   else if(prevDim < dim)
   {
      // adjust memory area to fit a new column, if needed
      int size = sizeof(double)*(dim-1);
      for(int row=dim-2; row >= 1; row--)
      {
         double *src  = Q + (row*(dim-1));
         double *dest = Q + (row*(dim-1))+ row;
         memmove(dest, src, size);
      }
   }
   
   prevDim = dim;
   
   // store new gradient and offset
   gradType = SML::DENSE;
   if(a.Density() < SPARSE_TO_DENSE_THRESHOLD)
      gradType = SML::SPARSE;
   
   if(gradientSet[idx])
      delete gradientSet[idx];
   gradientSet[idx] = new TheMatrix(a, gradType);
   
   offsetSet[idx] = b;
   
   // update time stamps
   enterTimeStamp[idx] = iter;
   activeTimeStamp[idx] = iter;  
   
   // update QP hessian matrix and vectors
   
   // Insert row into Q
   int theRow = idx*dim;
   for(int i=0; i < dim; i++) 
   {
      double value = 0;
      gradientSet[idx]->Dot(*gradientSet[i],value);
      Q[theRow + i] = value*QPScale/bmrmLambda;
   }
   
   // insert new column into Q[:,idx]
   for(int i=0; i < dim; i++)
   {
      Q[i*dim + idx] = Q[theRow + i];
   }
   
   if(useDiagQuadForm)
   {	
      // check dimensionality of B and feature vector
      if(a.Length() != B->Length())
      {
         ostringstream oss;
         oss << "Error: Feature vector (" << a.Length() << ") and B matrix (" << B->Length() << ") have incompatible dimensions" << endl;
         throw CBMRMException(oss.str(), "CL2N2_BMRMDualInnerSolver::Update()");
      }
      
      gradientSet[idx]->ElementWiseMult(*invB);
      double value = 0.0;
      gradientSet[idx]->Dot(a,value);
      Q[theRow+idx] = value*QPScale/bmrmLambda;
   }
   
   // update vectors      
   x[idx] = 0;
   f[idx] = -offsetSet[idx]*QPScale;
   
#ifdef QPDEBUG
   // Q matrix update correctness check!    
   MatrixCorrectnessCheck();
#endif  
}



int CL2N2_BMRMDualInnerSolver::AggregateGradients()
{
   // if aggregated gradient (AG) does not exist, let the earliest gradient be it
   if(aggGradIdx == -1)
   {
      aggGradIdx = 0;
      for(int i=0; i<dim; i++)
      {
         if(enterTimeStamp[aggGradIdx] > enterTimeStamp[i])
            aggGradIdx = i;
      }
        activeTimeStamp[aggGradIdx] = AGG_GRAD_TIME_STAMP;
   }
   // if AG had been moved to somewhere, find it
   else if(activeTimeStamp[aggGradIdx] != AGG_GRAD_TIME_STAMP)
   {
      for(int i=0; i<dim; i++)
         if(activeTimeStamp[i] == AGG_GRAD_TIME_STAMP)
         {
            aggGradIdx = i;
            break;
         }
   }
   // note that AG will almost never be deleted due to the large value of AGG_GRAD_TIME_STAMP
   
   
   // find most laziest/earliest gradient besides the agg. grad.
   int t1 = 0;  // for earliest gradient
   int t2 = 0;  // for laziest gradient
   int t3 = 0;  // for most active gradient
   
   t1 = t2 = t3 = (aggGradIdx+1)%dim;
   
   for(int i=0; i<dim; i++)
   {
      if(i == aggGradIdx) continue;
      t1 = (enterTimeStamp[t1] > enterTimeStamp[i] || t1 > i) ? i : t1;
      t2 = (activeTimeStamp[t2] > activeTimeStamp[i]) ? i : t2;
      t3 = (activeTimeStamp[t3] < activeTimeStamp[i]) ? i : t3;
   }
    
   // if laziest is not that lazy, we pick the earliest
   int idx = (t2 == t3) ? t1 : t2;
   
   assert(idx != aggGradIdx);
   if(idx == aggGradIdx) { cout << "idx == aggGradIdx" << endl; exit(EXIT_FAILURE); }
    
   // aggregation of gradients
   
   if(fabs(x[aggGradIdx]) < SML::ZERO_EPS)
   {
      // swap aggGradIdx with idx
      if(gradientSet[aggGradIdx]) 
         delete gradientSet[aggGradIdx];
      gradientSet[aggGradIdx] = gradientSet[idx];
      gradientSet[idx] = 0;
      
      offsetSet[aggGradIdx] = offsetSet[idx];
      offsetSet[idx] = 0;
      
      x[aggGradIdx] = x[idx];
   }
   else if(fabs(x[aggGradIdx]) > SML::ZERO_EPS && fabs(x[idx] > SML::ZERO_EPS))
   {
      gradientSet[aggGradIdx]->Scale(x[aggGradIdx]);      
      gradientSet[idx]->Scale(x[idx]);      
      gradientSet[aggGradIdx]->Add(*gradientSet[idx]);
      gradientSet[aggGradIdx]->Scale(1.0/(x[aggGradIdx]+x[idx]));
      
      offsetSet[aggGradIdx] *= x[aggGradIdx];
      offsetSet[aggGradIdx] += x[idx]*offsetSet[idx];
      offsetSet[aggGradIdx] /= (x[aggGradIdx]+x[idx]);
      
      x[aggGradIdx] = 0.0;
   }
   
   // update the hessian matrix and vectors
   for(int i=0; i<dim; i++)
   {
      // note that the entries Q[aggGradIdx][idx] and Q[idx][aggGradIdx] will be replaced correctly 
      //  when new gradient take the position #idx#
      if(i == idx) continue;
      
      double value = 0.0;
      gradientSet[aggGradIdx]->Dot(*gradientSet[i],value);
      Q[aggGradIdx*dim + i] = value*QPScale/bmrmLambda;
      Q[i*dim + aggGradIdx] = Q[aggGradIdx*dim + i];
    }

    f[aggGradIdx] = -offsetSet[aggGradIdx]*QPScale;     
    
    return idx;
}


/** Remove ALL idle gradients, if possible.
 *  Side effect:
 *  1.  change dim
 */
int CL2N2_BMRMDualInnerSolver::RemoveAllIdleGradients()
{
   int first_idle = 0;
   int last_idle = dim-1;
   int removeCnt = 0;
   
   // active gradients float to the top
   while(first_idle < last_idle)
   {                
      // look for the bottom most to-keep gradient and remove bottom most to-remove gradients
      while((last_idle > 0) && (iter-activeTimeStamp[last_idle] >= gradIdleAge)) 
      {
         if(gradientSet[last_idle]) 
         {
            delete gradientSet[last_idle];
            gradientSet[last_idle] = 0;                                
            offsetSet[last_idle] = 0;
            x[last_idle] = 0;
            f[last_idle] = 0; 
            activeTimeStamp[last_idle] = 0;
            enterTimeStamp[last_idle] = 0;
         }                        
         removeCnt++;
         last_idle--;
      }
      
      // look for top most to-remove gradient
      while((first_idle < dim) && (iter-activeTimeStamp[first_idle] < gradIdleAge)) 
         first_idle++;
      
      // replace the top most to-remove gradient with the bottom most to-keep gradient
      if(first_idle < last_idle)
      {       
         // 1. remove/replace the elements in gradientSet, offsetSet, x, f, activeTimeStamp                        
         if(gradientSet[first_idle]) 
            delete gradientSet[first_idle];
         gradientSet[first_idle] = gradientSet[last_idle];
         gradientSet[last_idle] = 0;
         
         offsetSet[first_idle] = offsetSet[last_idle];
         offsetSet[last_idle] = 0;
         
         x[first_idle] = x[last_idle];
         x[last_idle] = 0;
         
         f[first_idle] = f[last_idle];                                                     
         f[last_idle] = 0;
         
         activeTimeStamp[first_idle] = activeTimeStamp[last_idle];
         activeTimeStamp[last_idle] = 0;
         
         enterTimeStamp[first_idle] = enterTimeStamp[last_idle];
         enterTimeStamp[last_idle] = 0;
         
         // 2. memmove row                
         memmove(Q+first_idle*dim, Q+last_idle*dim, sizeof(double)*dim);
         
         // 3. memmove column
         for(int i=0; i<last_idle; i++)
            Q[i*dim + first_idle] = Q[i*dim + last_idle];                                
         
         // 4. one more to-remove gradient sink to the bottom
         removeCnt++;
      }
   }
   
   // make the hessian matrix compact
   if(removeCnt > 0)
   {
      assert(dim >= removeCnt);
      
      for(int i=1; i<dim; i++)
         memmove(Q+i*(dim-removeCnt), Q+i*dim, sizeof(double)*(dim-removeCnt));
      dim -= removeCnt;
   }
   dim++;
   
   return dim-1;  // i.e., use new slot
}


/** Remove the laziest gradient, if possible.
 *  Side effect:
 *  1. change dim
 */
int CL2N2_BMRMDualInnerSolver::RemoveLaziestIdleGradient()
{
   int idx = 0;
   
   for(int i=0; i < dim; i++) 
      idx = (activeTimeStamp[idx] > activeTimeStamp[i]) ? i : idx;
   
   // gradient not old enough
   if(iter - activeTimeStamp[idx] < gradIdleAge) 
   {
      idx = dim;
      dim++;
   }
   return idx;
}



/** Compute the value of regularizer
 *
 *  @param w [read] weight vector
 *  @return regularizer value
 */
double CL2N2_BMRMDualInnerSolver::ComputeRegularizerValue(TheMatrix &w)
{
   double regval = 0.0;
   if(useDiagQuadForm)
   {
      // this is inefficient; refine later
      TheMatrix Bw(w);
      Bw.ElementWiseMult(*B);
      Bw.Dot(w,regval);
   }
   else
   {
      w.Norm2(regval);
      regval *= regval;
   }
   regval *= bmrmLambda * 0.5; 
   return regval;
}


/** Return solution.
 *
 *  \param w      [write] solution vector (preallocated)
 *  \param objval [write] objective value
 */
void CL2N2_BMRMDualInnerSolver::GetSolution(TheMatrix& w, double &objval)
{
   assert(x != 0);
   double factor = 1.0/bmrmLambda;  // bmrmLambda is regularization constant
   
   double threshold = SML::ZERO_EPS;
   
   // compute objective value (explicitly)
   // n := dimensionality of x (the solution of QP)
   // Q := hessian matrix
   // f := linear part of obj
   
   objval = 0.0;
   double tmp = 0.0;
   double fx = 0.0;
   
   for(int i=0; i < dim; i++) 
   {
      tmp = 0;
      if(x[i] > threshold)
      {
         for(int j=0; j < dim; j++)
            if(x[j] > threshold)
               tmp += Q[j + i*dim]*x[j];
         objval += x[i]*tmp;
         fx += f[i]*x[i];
      }
   }        
   
   // since the dual of minimization problem is maximization
   // and we solved the negated maximization version,
   // so the objective value should be -objval
   objval = -0.5*objval - fx;
   objval /= QPScale;
   
   
   // Compute new w
   w.Zero();
   for(int i=0; i < dim; i++)
      if(x[i] > threshold)
         w.ScaleAdd(-x[i], *gradientSet[i]);
   w.Scale(factor);
   
   
   // update time-stamp of active gradients
   for(int i=0; i < dim; i++)
      if(x[i] > threshold && activeTimeStamp[i] <= iter)
         activeTimeStamp[i] = iter;         
   
   
   // count nonzero terms
   if(verbosity > 0) 
   {
      double smallest = 1e10;
      double sum_x = 0.0;
      int smallest_idx = -1;
      int howmany = 0;
      
      for(int i=0; i < dim; i++) 
      {
         if(x[i]> threshold) 
         {
            howmany++;
            if(smallest > x[i])
            {
               smallest = x[i];
               smallest_idx = i;
            }
            sum_x += x[i];
         }
      }
      cout << "qp size =  " << dim << endl;
      cout << "nnz(x) = " << howmany << endl;
      cout << "arg,min x = " << smallest_idx << ", " << smallest << endl << endl;
   }
}


void CL2N2_BMRMDualInnerSolver::Solve(TheMatrix& w, TheMatrix& grad, double loss, double &objval)
{
   double w_dot_grad = 0.0;
   w.Dot(grad, w_dot_grad);
   Update(grad, loss - w_dot_grad);
   
   SolveQP();
   GetSolution(w, objval);
   
   if(verbosity >= 3)
      for(int i=0; i<dim; i++)
         printf("x[%d] = %0.4e\n",i,x[i]);
}



#ifdef QPDEBUG
/** Check if the Q matrix update is correct
 */
void CL2N2_BMRMDualInnerSolver::MatrixCorrectnessCheck() 
{
   cout << "in Q matrix correctness check! " << endl;
   double *correctmat = new double[dim*dim];
   memset(correctmat, 0, sizeof(double)*dim*dim);
   
   for(int i=0; i < dim; i++)
   {
      for(int j=i; j < dim; j++)
      {
         double value = 0;
         gradientSet[i]->Dot(*gradientSet[j], value);
         correctmat[i*dim + j] = QPScale*value/bmrmLambda;
         correctmat[j*dim + i] = correctmat[i*dim+j];
      }
   }
   
   for(int i=0; i < dim*dim; i++)
   {
      if(fabs(correctmat[i] - Q[i]) > 1e-15) 
      {
         cout << "residual : " << fabs(correctmat[i] - Q[i]) << endl;
         cout << "Q matrix update correctness check : i = " << i << ", correct Q[i] = " << correctmat[i] << endl;
         cout << "ERROR: Q matrix is incorrect ! " << endl;
         exit(EXIT_FAILURE);
      }
   }
   
   if(correctmat) 
      delete [] correctmat;
}
#endif

#endif
