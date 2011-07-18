/* Copyright (c) 2009, INI
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
 * Authors: Jan Funke (funke@ini.phys.ethz.ch)
 *
 * Created: (07/14/2011)
 */

#include <iostream>

#include "common.hpp"
#include "softmarginloss.hpp"
#include "configuration.hpp"

using namespace std;

SoftMarginLoss::SoftMarginLoss(CModel* model, CConsVecData* data) :
	_data(data),
	_model(model) {

	cout << "[SoftMarginLoss] initialising..." << endl;

	_numFeatures  = _data->dim();
	_numVariables = _data->NumOfLabel();

	cout << "[SoftMarginLoss] dataset consists of " << _numFeatures
	     << " features and " << _numVariables << " variables." << endl;

	// initialize model
	if(!_model->IsInitialized())
		_model->Initialize(1, _numFeatures, _data->bias());

}

void
SoftMarginLoss::ComputeLoss(double& loss) {

	cout << "[SoftMarginLoss::ComputeLoss] ..." << endl;

	loss = 0;
}

void
SoftMarginLoss::ComputeLossAndGradient(double& loss, TheMatrix& grad) {

	cout << "[SoftMarginLoss::ComputeLossAndGradient] ..." << endl;

	TheMatrix& weights = _model->GetW();
	TheMatrix f(_numVariables, 1);
	_data->XMultW(weights, f);

	f.Print();

	loss = 0;
}
