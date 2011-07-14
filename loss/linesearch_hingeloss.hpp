/* Copyright (c) 20069, NICTA
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
 * Created: 16/05/2008
 *
 * Last Updated: 16/01/2009
 */

#ifndef _LINESEARCH_HINGELOSS_HPP_
#define _LINESEARCH_HINGELOSS_HPP_

#include <vector>
#include "sml.hpp"
#include "binaryclassificationloss.hpp"
#include "model.hpp"
#include "hinge_linesearch.hpp"


/**  Class to represent binary hinge loss function with linesearch procedure:
 *   loss = max(0, margin - y*f(x))
 *   where f(x) := <w, x> and margin \geq 0. By default margin = 1.
 *
 *   Note that the loss definition is exactly the same as hingeloss but 
 *   we exploit some precomputed values here which 
 */
class CLinesearch_HingeLoss : public CBinaryClassificationLoss 
{
   protected:
      /** For line-search 
       */
      CHinge_Linesearch *ls;  
      TheMatrix *wb;
      TheMatrix *p;
      TheMatrix *fb;
      TheMatrix *df;
      TheMatrix *delta;
      double bmrmLambda;
      double theta; // adjustment between w_{t+1}^b and w_t
      double loss_of_wb;
      double eta;   // line search step size
      
      // not used in this line search based loss; but still need to be declared
      void Loss(double& loss, TheMatrix& f) {}
      void LossAndGrad(double& loss, TheMatrix& f, TheMatrix& l) {}
      
   public:    
      
      CLinesearch_HingeLoss(CModel* &model, CVecData* &data);
      virtual ~CLinesearch_HingeLoss();
      virtual void ComputeLoss(double &loss);
      virtual void ComputeLossAndGradient(double &loss, TheMatrix &grad);
      virtual double GetLossOfWbest()
      {
         return loss_of_wb;
      }
      
      virtual void GetWbest(TheMatrix &w_return)
      {
         w_return.Assign(*wb);
      }
      virtual double GetLineSearchStepSize()
      {
         return eta;
      }
};

#endif
