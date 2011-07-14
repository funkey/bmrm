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
 * Last Updated: (19/11/2007)   
 */

#ifndef _UNIVARIATEREGRESSIONLOSS_HPP_
#define _UNIVARIATEREGRESSIONLOSS_HPP_

#include "sml.hpp"
#include "scalarloss.hpp"
#include "model.hpp"

/**  
 * Abstract base class for encapsulating univariate regression losses. 
 * 
 * The assumption is that the labels are real valued and the loss measures
 * the discrepancy between predicted values and the real valued
 * labels. This class also contains prediction and performance
 * evaluation metrics for univariate regression loss. 
 */
class CUnivariateRegressionLoss : public CScalarLoss 
{
   public:    
      CUnivariateRegressionLoss(CModel* &model, CVecData* &data)
         : CScalarLoss(model, data) {}
      virtual ~CUnivariateRegressionLoss() {}
      
   protected:
      virtual void Loss(double& loss, TheMatrix& f) = 0;
      virtual void LossAndGrad(double& loss, TheMatrix& f, TheMatrix& l) = 0;
      void Transform(TheMatrix& f);
      void Perf(TheMatrix& f, TheMatrix& predict);
};

#endif
