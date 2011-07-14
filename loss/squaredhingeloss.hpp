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
 * Created: (02/11/2007) 
 *
 * Last Updated: (12/01/2009)   
 */

#ifndef _SQUAREDHINGELOSS_HPP_
#define _SQUAREDHINGELOSS_HPP_

#include "sml.hpp"
#include "binaryclassificationloss.hpp"
#include "model.hpp"

/**  
 * Class to represent binary squared hinge loss function:
 * loss = 0.5*max(0, margin - y*f(x))^2
 * where f(x) := <w, x> and margin \geq 0. By default margin = 1.
 */
class CSquaredHingeLoss : public CBinaryClassificationLoss 
{
   protected:
      void Loss(double& loss, TheMatrix& f);
      void LossAndGrad(double& loss, TheMatrix& f, TheMatrix& l);
      
   public:    
      CSquaredHingeLoss(CModel* &model, CVecData* &data)
         : CBinaryClassificationLoss(model, data) {}
      virtual ~CSquaredHingeLoss() {}
};

#endif
