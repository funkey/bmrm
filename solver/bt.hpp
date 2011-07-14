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

#ifndef _BT_HPP_
#define _BT_HPP_

#include "common.hpp"
#include "sml.hpp"   // the matrix library wrapper
#include "model.hpp"
#include "loss.hpp"
#include "regularizer.hpp"
#include "solver.hpp"

/**   Class for BT solver.
 */
class CBT : public CSolver
{
   public:
      CBT(CModel *model, CLoss *loss, CRegularizer *reg);
      CBT(CModel *model, CLoss *loss);
      virtual ~CBT();
      virtual void Train();
      
   protected:
      void ReadParams();
      void ModelLossAndGrad(const int& n,
                            const double* w,
                            double& obj,
                            double* grad);
      
      // // model
      // CModel *_model;
      
      // // loss function
      // CLoss *_loss;

      // regularizer
      CRegularizer *_reg;

      // gradient vector for loss function
      TheMatrix *my_grad;
      
      // maximum iterations 
      int max_iter;      
      
      // max number of function/gradient evaluations
      int max_eval;  
      
      // regularization constant
      double lambda;            
      
      // debug level
      int verb;
      
      // buffer size
      int buf_size;     
      
      // error tolerance
      double tol;  
      
      // lower bound on objective
      double fm;    
};

#endif
