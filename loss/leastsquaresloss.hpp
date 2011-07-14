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
 * Last Updated: 15/01/2009
 */

#ifndef _LEASTSQUARESLOSS_HPP_
#define _LEASTSQUARESLOSS_HPP_

#include "sml.hpp"
#include "univariateregressionloss.hpp"
#include "model.hpp"

/**  
 * Class to represent least squares regression loss function:
 * loss = (y-<w,x>)^2
 */
class CLeastSquaresLoss : public CUnivariateRegressionLoss
{
   protected:
      void Loss(double& loss, TheMatrix& f);
      void LossAndGrad(double& loss, TheMatrix& f, TheMatrix& l);
      
   public:    
      CLeastSquaresLoss(CModel* &model, CVecData* &data)
         : CUnivariateRegressionLoss(model, data) {}
      virtual ~CLeastSquaresLoss() {}
};

#endif
