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
 * Created      : 25/07/2007
 * Last Updated :
 */

/* 
 * Purpose      : L1-norm Linear Program solver using COIN-OR Clp 
 *
 */

#ifndef _L1N1_CLP_CPP_
#define _L1N1_CLP_CPP_

#include "l1n1_clp.hpp"
#include "configuration.hpp"
#include "bmrmexception.hpp"

// Clp headers
#include "CoinBuild.hpp"


using namespace std;

/*
    L1-Norm linear program:

    minimize   :  d'x
       x
     
    subject to : | -1_cm  A_mn  -A_mn | |\xi|    
                                        | u | <= |-B| 
                                        | v |    

		 x >= 0

    where      : d := [1; lambda; -lambda] (length: [1,n,n])
                 C is a positive scalar
                 lambda is a column vector of size n
                 x := [\xi; u; v]
                 u,v column vectors of size n
                 1_c# is a column vector of 1 of size #,
                 A_#* is a dense matrix of size #-by-*, and
                 m << n (m starts from 1 and increased by 1 after 
                 each lp optimization round.)
		 B is a column vector of size m
		 w = u-v
*/


CL1N1_Clp::CL1N1_Clp(double lambda, const int &thedim)
   : CBMRMInnerSolver(lambda),
     LPScale(1.0),
     newRowIndices(0),
     newRowElements(0),
     gradIdleAge(10),
     removeAllIdleGrad(false)
{
   // parameters
   iter = 0;
   dim = thedim;
   numOfConstraint = 0;
   
   // sanity check
   assert(dim > 0);
   
   // read configuration
   Configuration &config = Configuration::GetInstance();
   if(config.IsSet("L1N1_Clp.gradIdleAge"))
   {
      gradIdleAge = config.GetInt("L1N1_Clp.gradIdleAge");
      gradIdleAge = max(2,gradIdleAge);
   }
   
   if(config.IsSet("L1N1_Clp.removeAllIdleGradients"))
      removeAllIdleGrad = config.GetBool("L1N1_Clp.removeAllIdleGradients");
   
   if(config.IsSet("L1N1_Clp.scale"))
      LPScale = config.GetDouble("L1N1_Clp.scale");
   
   // build simplex model
   sim = new ClpSimplex();
   sim->setLogLevel(0);
   sim->resize(0, 1+2*dim);
   
   // set objective coefficients    
//   sim->setObjCoeff(0, 1.0);
//   for(int i=1; i<2*dim+1; i++) 
//       sim->setObjCoeff(i, bmrmLambda);
   
   sim->setObjCoeff(0, LPScale/bmrmLambda);
   for(int i=1; i<2*dim+1; i++) 
      sim->setObjCoeff(i, LPScale);
   
   // set bounds for variables (i.e., columns as in COIN Clp terminology)
   if(nonNegativeSlack)
      sim->setColumnBounds(0, 0.0, COIN_DBL_MAX);
   else
      sim->setColumnBounds(0, -COIN_DBL_MAX, COIN_DBL_MAX);  
   
   for (int i=1; i<2*dim+1; i++)
      sim->setColumnBounds(i, 0.0, COIN_DBL_MAX);
   
   
   // something useful in constraint matrix update phase
   newRowIndices = (int*)calloc(2*dim+1, sizeof(int));
   newRowElements = (double*)calloc(2*dim+1, sizeof(double));    
   
   for(int i=0; i<2*dim+1; i++)
      newRowIndices[i] = i;
   
   newRowElements[0] = -1.0;
   // Update() set the rest of the elements in newRowElements[]
   
}


CL1N1_Clp::~CL1N1_Clp()
{        
   if(sim) delete sim;
   if(newRowIndices) free(newRowIndices);
   if(newRowElements) free(newRowElements);  
}


void CL1N1_Clp::Solve(TheMatrix& w, TheMatrix& grad, double loss, double &objval)
{
   double w_dot_grad = 0.0;
   w.Dot(grad, w_dot_grad);
   Update(grad, loss - w_dot_grad);
   
   // update lower bound on \xi
   double *solution = sim->primalColumnSolution(); 
   double lowerbound = solution[0];
   sim->setColumnBounds(0, lowerbound, COIN_DBL_MAX);
   
   //sim->dual();
   sim->primal();
   
   GetSolution(w, objval);
}


/** Clear constraint matrix
 */
void CL1N1_Clp::Reset()
{
   // delete rows in constraint matrix
   int nrows = sim->numberRows();
   int *rows = new int[nrows];
   for(int i=0; i<nrows; i++) 
      rows[i] = i;
   sim->deleteRows(nrows, rows);      
   delete rows;
}

void CL1N1_Clp::ChangeLambda(double newLambda)
{
   bmrmLambda = newLambda;
   sim->setObjCoeff(0, LPScale/bmrmLambda);
}


/**  Update matrix and vectors needed in optimisation round.
 *
 *   \param a [read] Gradient
 *   \param b [read] Offset 
 */
void CL1N1_Clp::Update(TheMatrix& a, double b)
{
   if(verbosity > 1)
   {
      double anorm = 0.0;
      a.Norm2(anorm);
      cout << "|a|_2 = " << anorm << endl;
   }
   
   iter++;
   
   int length = a.Length();
   assert(length == dim);
   
   // reset the content of newRowElements and newRowIndices
   memset(newRowElements, 0, sizeof(double)*(1+dim+dim));
   memset(newRowIndices, 0, sizeof(int)*(1+dim+dim));
   newRowElements[0] = -1;
   newRowIndices[0] = 0;
   
   double val = 0;
   unsigned int nnz = 0;
   
   for(int i=0; i<length; i++)  
   {
      a.Get(i,val);
      if(fabs(val) > SML::ZERO_EPS) 
      {
         newRowElements[1+nnz] = val;
         newRowIndices[1+nnz] = i+1;  // gradient starts from second column
         nnz++;
      }
   }
   
   for(unsigned int i=0; i<nnz; i++)  
   {
      newRowElements[1+nnz+i] = -newRowElements[1+i];
      newRowIndices[1+nnz+i]  = newRowIndices[1+i]+dim;
   }
   
   sim->addRow(1+nnz+nnz, newRowIndices, newRowElements, -COIN_DBL_MAX, -b);
   numOfConstraint++;
   activeTimeStamp.push_back(iter);
   
   if(verbosity > 0) 
   {
      cout << "num cols: " << sim->numberColumns() << endl;
      cout << "num rows: " << sim->numberRows() << endl;
   }
}


double CL1N1_Clp::ComputeRegularizerValue(TheMatrix &w)
{
    double regval = 0.0;
    w.Norm1(regval);
    regval *= bmrmLambda;
    return regval;
}


/** Return solution.
 *
 *  \param w      [write] solution vector (preallocated)
 *  \param objval [write] objective value
 */
void CL1N1_Clp::GetSolution(TheMatrix& w, double &objval)
{
   double *solution = 0;
   solution = sim->primalColumnSolution(); 
   assert(solution);
   
   // construct new w
   unsigned int length = w.Length();
   for(unsigned int i=0; i<length; i++) 
      w.Set(i, solution[1+i] - solution[1+dim+i]);   // w := u-v
   
   //objval = sim->objectiveValue();
   objval = bmrmLambda/LPScale*sim->objectiveValue();
   
   if(verbosity > 1)
   {
      cout << "LP obj val = " << objval << endl;
      cout << "xi = " << solution[0] << endl;
   }
   
   vector<int> del;
   int nrows = sim->numberRows();
   
   if(removeAllIdleGrad)
   {
      for(int i=0; i<nrows; i++)
      {
         if(sim->getRowStatus(i)==ClpSimplex::atUpperBound)
            activeTimeStamp[i] = iter;
         
         if(iter-activeTimeStamp[i] >= gradIdleAge)
            del.push_back(i);
      }
      
      if(verbosity > 0)
         cout << "num of deleted constraints: " << del.size() << endl;
      
      sim->deleteRows(del.size(), &(del[0]));
      int size=del.size();
      for(int i=size-1; i>0; i--)
         activeTimeStamp.erase(activeTimeStamp.begin() + del[i]);
      numOfConstraint -= del.size();
   }
   else
   {
      int oldestTS = activeTimeStamp[0];
      int oldestIdx = 0;
      int nrows = sim->numberRows();
      
      for(int i=0; i<nrows; i++)
      {
         if(sim->getRowStatus(i)==ClpSimplex::atUpperBound)
            activeTimeStamp[i] = iter;
         
         if(activeTimeStamp[i] < oldestTS)
         {
            oldestTS = activeTimeStamp[i];
            oldestIdx = i;
         }
      }
      
      // delete oldest inactive constraint
      if(iter-oldestTS >= gradIdleAge)
      {
         if(verbosity > 0)
            cout << "deleted constraint: " << oldestIdx << endl;
         
         sim->deleteRows(1,&oldestIdx);
         activeTimeStamp.erase(activeTimeStamp.begin() + oldestIdx);
         numOfConstraint--;
      }
   }
}

#endif
