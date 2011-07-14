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
 * Created: 02/11/2007
 *
 * Last Updated:
 */

#ifndef _EPSILONINSENSITIVELOSS_HPP_
#define _EPSILONINSENSITIVELOSS_HPP_

#include "sml.hpp"
#include "univariateregressionloss.hpp"
#include "model.hpp"

/**  
 * Class to represent epsilon insensitive regression loss:
 * loss = max(0, |y - f(x)| - epsilon)
 * where f(x) := <w, x> and epsilon \geq 0. By default margin = 0.1.
 */
class CEpsilonInsensitiveLoss : public CUnivariateRegressionLoss 
{
   protected:
      /** The epsilon parameter in the \epsilon-insensitve loss function
       */
      double epsilon;
      
      void Loss(double& loss, TheMatrix& f);
      void LossAndGrad(double& loss, TheMatrix& f, TheMatrix& l);
      
   public:    
      CEpsilonInsensitiveLoss(CModel* &model, CVecData* &data);
      virtual ~CEpsilonInsensitiveLoss() {}
};

#endif
