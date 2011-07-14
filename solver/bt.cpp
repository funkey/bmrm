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
 * Authors: S V N Vishwanathan (vishy@stat.purdue.edu)
 *          Nic Schraudolph (nic@schraudolph.org)
 * Created: (19/11/2008) 
 *
 */

#ifndef _BT_CPP_
#define _BT_CPP_

#include "common.hpp"
#include "timer.hpp"
#include "info.hpp"
#include "configuration.hpp"
#include "regularizerfactory.hpp"
#include "bt.hpp"

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

using namespace std;

// keep your hands off this. This is only for use by the loss and grad function. 
CLoss* _my_loss;
CRegularizer* _my_reg;
CModel* _my_model;
TheMatrix* _my_grad;
double _my_lambda;
  
extern "C"
{
  void bt_(const int&, double[], double&, const double[],
           void(const int&, const double*, double&, double*),
           const double&, const double&, const int&, const int&, const int&,
           const int&, int&, const int[], double[], const int&);
}

CBT::CBT(CModel *model, CLoss *loss, CRegularizer *reg)
   :CSolver(model,loss),
    _reg(reg)
{ 
   my_grad = new TheMatrix(model->GetW());
   _my_grad = my_grad;
   _my_loss = loss;
   _my_reg = reg;
   _my_model = model;
}

CBT::CBT(CModel *model, CLoss *loss)
   : CSolver(model,loss)
{ 
   _reg = CRegularizerFactory::GetRegularizer();
   my_grad = new TheMatrix(model->GetW());
   _my_grad = my_grad;
   _my_loss = loss;
   _my_reg = _reg;
   _my_model = model;
}

CBT::~CBT()
{
   if(my_grad) delete my_grad;
   if(_reg) delete _reg;
}


// Calculates reg+loss and its gradient
// Very very very fragile and hacky

void bt_func_and_grad(const int& n,
                      const double* x,
                      double& obj,
                      double* grad)
{
   obj = 0.0;
   double val = 0.0;
   double gi = 0.0;

   // update w in model with new solution x
   TheMatrix *my_w = &_my_model->GetW();
   for (int i = 0; i < n; i++)
      my_w->Set(i, x[i]);

   _my_reg->ComputeRegAndGradient(*_my_model, val, *_my_grad);
   for(int i=0; i<n; i++)
   {
      _my_grad->Get(i,gi);
      grad[i] = _my_lambda*gi;
   }
   obj = _my_lambda*val;

   _my_loss->ComputeLossAndGradient(val, *_my_grad);
   for(int i=0; i<n; i++)
   {
      _my_grad->Get(i,gi);
      grad[i] += gi;
   }
   obj += val;
}


void CBT::Train()
{
   TheMatrix &w = _model->GetW();

   // Read in optimizer parameters
   ReadParams();
  
   // svnvish: BUGBUG
   // internals - do not edit below here
   int *iwork = new int[buf_size];
   const int n = w.Length(); 
   const int lwork = (buf_size + n + 9)*buf_size + 4*n + 10;
   double *work = new double[lwork];   
   int ifail = 0;            // error code
   double f = 0;                 // function value
   double *g = new double[n];              // gradient
   double *x = new double[n];              // parameter (temporary)
   
   for (int i = 0; i < n; i++) 
   {
      double tmp = 0;
      w.Get(i, tmp);
      x[i] = 0;
   }
  
   bt_func_and_grad(n, x, f, g);   // to calculate intial f and g
   
   bt_(n, x, f, g, bt_func_and_grad, fm, tol, 2*max_iter, max_iter,
       buf_size, verb, ifail, iwork, work, lwork); 

   if(iwork) delete iwork;
   if(work) delete work;
   if(g) delete g;
   if(x) delete x;
   
   if (ifail != 0) // cerr << isinf::argv0 << "bt: ";
      cerr << "bt: ";
   switch (ifail){
      case 5:
         cerr << "QP solver failed" << endl;
         exit(ifail);
      case 3:
         cerr << "fm too large" << endl;
         exit(ifail);
      case 2:
         cerr << "warning: maxit exceeded" << endl;
         break;
      case 1:
         cerr << "warning: maxcom exceeded" << endl;
         break;
      case 0:
         break;
      default:
         cerr << "unknown error " << ifail << endl;
         exit(ifail);
   }
   
   // return f;  
}

 

// Read and store relevant parameters for the optimizer
void CBT::ReadParams()
{
   // make sure configuration file is read before this!
  Configuration &config = Configuration::GetInstance();  
  
  if(config.IsSet("BT.maxNumOfIter")){
     max_iter = config.GetInt("BT.maxNumOfIter");
     if(max_iter < 0)
        throw CBMRMException("BT.maxNumOfIter must be > 0\n",
                             "CBT::ConfirmProgramParameters()");
  }else{
     max_iter = 1000;
  }
  
  //svnvish: BUGBUG
  // arbitrary
  max_eval = 2*max_iter;
  
  if(config.IsSet("BMRM.lambda")){
    lambda = config.GetDouble("BMRM.lambda");
    if(lambda <= 0)
       throw CBMRMException("BMRM.lambda must be > 0\n",
                            "CBT:ConfirmProgramParameters()");
  }else{
     throw CBMRMException("BMRM.lambda must be specified \n",
                          "CBT:ConfirmProgramParameters()");
  }
  
  // svnvish: BUGBUG
  // hack
  _my_lambda = lambda;

  if(config.IsSet("BT.debugLevel")){
    verb = config.GetInt("BT.debugLevel");
    if(verb < 0)
       throw CBMRMException("BT.debugLevel must be >= 0\n",
                            "CBT:ConfirmProgramParameters()");
  }else{
     verb = 2;
  }
  
  if(config.IsSet("BT.bufferSize")){
     buf_size = config.GetInt("BT.bufferSize");
     if(buf_size <= 0)
        throw CBMRMException("BT.bufferSize must be > 0\n",
                           "CBT:ConfirmProgramParameters()");
  }else{
     buf_size = 10;
  }
  
  if(config.IsSet("BT.epsilonTol")){
     tol = config.GetDouble("BT.epsilonTol");
     if(tol <= 0)
        throw CBMRMException("BT.epsilonTol must be > 0\n",
                             "CBT:ConfirmProgramParameters()");
  }else{
     tol = 1e-6;
  }

  fm = 0.0;  
  if(config.IsSet("BT.fm"))
  {
     fm = config.GetDouble("BT.fm");
  }
  
}

#endif

