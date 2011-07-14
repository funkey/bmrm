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
 * Created: 18/02/2009
 *
 * Last Updated:
 */

#ifndef _HUBERHINGELOSS_HPP_
#define _HUBERHINGELOSS_HPP_

#include "sml.hpp"
#include "binaryclassificationloss.hpp"
#include "model.hpp"

/**  
 * Class to represent Huber hinge loss:
 *            / 0                       if yf > 1+h
 *  loss_i = <  ((1+h-yf)^2)/(4*h)  if |1-yf| <= h
 *            \ 1 - yf                  if yf < 1-h
 *
 * where f := <w, x>, h is user defined and typically in [0.01,0.5]
 * 
 * Reference: 
 * O. Chaplle 
 * Training a Support Vector Machine in the Primal
 * Neural Computation, 2007.
 */
class CHuberHingeLoss : public CBinaryClassificationLoss
{
   protected:
      /** The huber parameter
       */
      double h;
      
      void Loss(double& loss, TheMatrix& f);
      void LossAndGrad(double& loss, TheMatrix& f, TheMatrix& l);
      
   public:    
      CHuberHingeLoss(CModel* &model, CVecData* &data);
      virtual ~CHuberHingeLoss() {}
};

#endif
