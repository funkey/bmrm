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
	_model(model),
	_numFeatures(_data->dim()),
	_numVariables(_data->NumOfLabel()),
	_numEqConstraints(_data->NumOfEqualities()),
	_numIneqConstraints(_data->NumOfInequalities()),
	_solver(_numVariables, _numEqConstraints, _numIneqConstraints),
	_costFactor(1.0),
	_linearCostContribution(_numVariables, 1, SML::DENSE) {

	cout << "[SoftMarginLoss] initialising..." << endl;

	cout << "[SoftMarginLoss] dataset consists of " << _numFeatures
	     << " features and " << _numVariables << " variables." << endl;

	// initialize model
	if(!_model->IsInitialized())
		_model->Initialize(1, _numFeatures, _data->bias());

	// set linear constraints
	_solver.SetEqualities(
			_data->GetEqualityCoefs(),
			_data->GetEqualityValues());
	_solver.SetInequalities(
			_data->GetInequalityCoefs(),
			_data->GetInequalityValues());

	// calculate cost contribution (does not change anymore)
	ComputeCostContribution(_data->labels());
}

void
SoftMarginLoss::ComputeLoss(double& loss) {

	cout << "[SoftMarginLoss::ComputeLoss] ..." << endl;

	TheMatrix grad(_numFeatures, 1, SML::DENSE);

	ComputeLossAndGradient(loss, grad);
}

void
SoftMarginLoss::ComputeLossAndGradient(double& loss, TheMatrix& grad) {

	cout << "[SoftMarginLoss::ComputeLossAndGradient] ..." << endl;

	// compute coefficients of the objective
	TheMatrix& weights = _model->GetW();
	TheMatrix f(_numVariables, 1);
	_data->XMultW(weights, f);

	// set objective in linear solver
	_solver.SetObjective(f, 1.0);

	// solve the ILP

	// get value of the objective
	loss = 0;

	// compute gradient
}

void
SoftMarginLoss::ComputeCostContribution(const TheMatrix& groundTruth) {

	cout << "[SoftMarginLoss::ComputeCostContribution] "
	     << "computing loss contributions for " << endl;
	groundTruth.Print();

	_constantCostContribution = 0;

	for (int i = 0; i < _numVariables; i++) {

		double value;

		// accumulate constant cost term
		groundTruth.Get(i, value);
		_constantCostContribution += value*_costFactor;

		// calculate linear cost coefficients
		_linearCostContribution.Set(i, (1.0 - 2.0*value)*_costFactor);
	}

	cout << "[SoftMarginLoss::ComputeCostContribution] "
	     << "constant contribution is " << _constantCostContribution
	     << endl;

	cout << "[SoftMarginLoss::ComputeCostContribution] "
	     << "linear contribution is " << endl;
	_linearCostContribution.Print();
}
