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

#ifndef _LEASTABSDEVLOSS_HPP_
#define _LEASTABSDEVLOSS_HPP_

#include "sml.hpp"
#include "univariateregressionloss.hpp"
#include "model.hpp"

/**  
 * Class to represent least absolute deviation regression loss function:
 * loss = abs(y-<w,x>)
 */
class CLeastAbsDevLoss : public CUnivariateRegressionLoss
{
   protected:
      void Loss(double& loss, TheMatrix& f);
      void LossAndGrad(double& loss, TheMatrix& f, TheMatrix& l);
      
   public:    
      CLeastAbsDevLoss(CModel* &model, CVecData* &data)
         : CUnivariateRegressionLoss(model, data) {}
      virtual ~CLeastAbsDevLoss() {}
};

#endif
