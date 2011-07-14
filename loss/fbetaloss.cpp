/* Copyright (c) 2009, NICTA
 * All rights reserved. 
 * 
 * The contents of this file are subject to the Mozilla Public License 
 * Version 1.1 (the "License"); you may not use this file except in 
 * compliance with the License. You may obtain a copy of the License at 
 * http://www.mozilla.org/MPL/ 
 * 
 * Software distributed under the License is distributed on an "AS IS" 
 * basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the 
 * License for the specific language governing rights and limitations 
 * under the License. 
 * 
 * Authors: Choon Hui Teo (ChoonHui.Teo@anu.edu.au)
 *
 * Created: (12/01/2009) 
 *
 * Last Updated:
 */

#ifndef _FBETALOSS_CPP_
#define _FBETALOSS_CPP_

#include "fbetaloss.hpp"
#include "sml.hpp"
#include "configuration.hpp"

using namespace std;


CFBetaLoss::CFBetaLoss(CModel* &model, CVecData* &data)
   : CBinaryClassificationLoss(model, data),
     m_pos(0),
     m_neg(0),
     orig_pidx(0),
     orig_nidx(0),
     fpos(0),
     fneg(0),
     beta(1.0)
{
   const TheMatrix &Y = _data->labels();
   double y_i = 0.0;
   for(unsigned int i=0; i < _data->size(); i++)
   {
      Y.Get(i,y_i);
      if(y_i > 0.0) 
         m_pos++;
      else
         m_neg++;
   }
   assert(m == (m_pos + m_neg));
   
   orig_pidx = new int[m_pos];        
   orig_nidx = new int[m_neg];        
   int p = 0;
   int n = 0;
   for(unsigned int i=0; i<m; i++) 
   {
      Y.Get(i,y_i);
      if(y_i > 0.0)
         orig_pidx[p++] = i;
      else
         orig_nidx[n++] = i;
   }
   pidx = new int[m_pos];
   nidx = new int[m_neg];
   fpos = new double[m_pos+1];
   fneg = new double[m_neg+1];

   // when there is no special scaling applied to the labels 
   // (check svm_struct_api.c::read_struct_examples()),
   // the relation between the C parameter of svmperf F1 score and
   // the lambda in bmrm F1 score is: lambda = 1/C
   scalingFactor = 100.0/m;
   
   Configuration &config = Configuration::GetInstance();
   if(config.IsSet("F_BETA.beta"))
   {
      beta = config.GetDouble("F_BETA.beta");
      if(beta <= 0.0)
         throw CBMRMException("Error: F_BETA.beta must be > 0.0","CFBetaLoss::CFBetaLoss()");
   }

   if(verbosity > 0)
   {
      cout << "In CFBetaLoss" << endl;
      cout << "Num. of +ve labels: " << m_pos << endl;
      cout << "Num. of -ve labels: " << m_neg << endl;
      cout << "Beta: " << beta << endl;
   }
}
        
        
CFBetaLoss::~CFBetaLoss()
{
   if(orig_pidx) delete [] orig_pidx;
   if(orig_nidx) delete [] orig_nidx;
   if(pidx) delete [] pidx;
   if(nidx) delete [] nidx;
   if(fpos) delete [] fpos;
   if(fneg) delete [] fneg;
}


/**  
 *  Compute FBeta loss. CAUTION: f is passed by reference and is
 *  changed within this function. This is done for efficiency reasons,
 *  otherwise we would have had to create a new copy of f. 
 *   
 *  @param loss [write] loss value computed.
 *  @param f [read/write] prediction vector. 
 */
void CFBetaLoss::Loss(double& loss, TheMatrix& f)
{
   LossAndGrad(loss, f, *l);
}


/**  
 *  Compute loss and partial derivative of FBeta loss w.r.t f
 *   
 *  @param loss [write] loss value computed.
 *  @param f [r/w] = X*w
 *  @param l [write] partial derivative of loss w.r.t. f
 */
void CFBetaLoss::LossAndGrad(double& loss, TheMatrix& f, TheMatrix& l)
{
   const TheMatrix &Y = _data->labels();
   double *f_array = f.Data();

   indirect_greater_than<double> igt(f_array); // in descending order   
   indirect_less_than<double> ilt(f_array); // in ascending order

   memcpy(pidx, orig_pidx, sizeof(int)*m_pos);
   sort(pidx, pidx+m_pos, igt);

   memcpy(nidx, orig_nidx, sizeof(int)*m_neg);
   sort(nidx, nidx+m_neg, ilt);

   fpos[0] = 0.0;
   for(unsigned int i=0; i<m_pos; i++)
      fpos[i+1] = fpos[i] + 2*f_array[pidx[i]];
   
   fneg[0] = 0.0;
   for(unsigned int i=0; i<m_neg; i++)
      fneg[i+1] = fneg[i] + 2*f_array[nidx[i]];

   loss = 0.0;
   l.Assign(Y);
   l.Scale(-2.0);
   double tmploss = 0;
   unsigned int best_tp = 0;
   unsigned int best_tn = 0;

   for(unsigned int tp=0; tp <= m_pos; tp++)
   {
      for(unsigned int tn=0; tn <= m_neg; tn++)
      {
         tmploss = Delta(tp,tn) - (fpos[m_pos]-fpos[tp]) + (fneg[m_neg]-fneg[tn]);
         if(tmploss > loss)
         {
            loss = tmploss;
            best_tp = tp;
            best_tn = tn;
         }
      }
   }

   for(unsigned int i=0; i < best_tp; i++)
      l.Set(pidx[i],0);
   
   for(unsigned int i=0; i < best_tn; i++)
      l.Set(nidx[i],0);
}

#endif
