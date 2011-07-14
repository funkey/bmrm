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
 * Created: (21/10/2008) 
 *
 * Last Updated:
 */

#ifndef _LOGISTICLOSS_HPP_
#define _LOGISTICLOSS_HPP_

#include <vector>
#include "sml.hpp"
#include "binaryclassificationloss.hpp"
#include "model.hpp"


/** Class for Logistic loss function
 *  loss_i = ln(1+exp(-y_i*f_i))
 *  where f_i := <w,x_i>
 *
 *  NOTE: to address the overflow issue:
 *            /  ln(2)               if yf = 0  
 *  loss_i = <   ln(1+exp(-yf))      if yf > 0
 *            \  ln(1+exp(yf)) + yf  if yf < 0
 */

class CLogisticLoss : public CBinaryClassificationLoss 
{
protected:
	void Loss(double& loss, TheMatrix& f);
	void LossAndGrad(double& loss, TheMatrix& f, TheMatrix& l);
	
public:    
	
	CLogisticLoss(CModel* &model, CVecData* &data)
		: CBinaryClassificationLoss(model, data) {}
	virtual ~CLogisticLoss() {}
};

#endif
