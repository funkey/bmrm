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
	if (_numEqConstraints > 0)
		_solver.SetEqualities(
				_data->GetEqualityCoefs(),
				_data->GetEqualityValues());
	if (_numIneqConstraints > 0)
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


	/////////////////////////////////////////////////////////////////
	// compute coefficients f and constant part c of the objective //
	/////////////////////////////////////////////////////////////////

	// the current weights
	TheMatrix& w = _model->GetW();

	cout << "[SoftMarginLoss::ComputeLossAndGradient] "
	     << "current weights: ";
	w.Print();

	// f = Xw (intermediate result)
	TheMatrix f(_numVariables, 1);
	_data->XMultW(w, f);

	// c = c1 - c2 = <y',1*gamma> - <y',Xw>

	// c1 = <y',1*gamma>
	double c1 = _constantCostContribution;

	// c2 = <y',Xw>
	double c2;
	_data->labels().Dot(f, c2);

	double c = c1 - c2;

	// f = Xw + (1-2y')*gamma
	f.Add(_linearCostContribution);

	// set objective in linear solver
	_solver.SetObjective(f, c, LinearProgramSolver::MAXIMIZE);

	///////////////////
	// solve the ILP //
	///////////////////

	// the solution vector
	TheMatrix y(_numVariables, 1, SML::SPARSE);

	// the return message of the solver
	string msg;

	bool success = _solver.Solve(y, loss, msg);

	cout << "[SoftMarginLoss::ComputeLossAndGradient] "
	     << "most offending solution: ";
	y.Print();
	cout << "[SoftMarginLoss::ComputeLossAndGradient] "
	     << "loss: " << loss << endl;

	if (!success) {

		cout << "[SoftMarginLoss::ComputeLossAndGradient]" << endl
		     << " **************** WARNING *******************" << endl
		     << " **** linear solver did not find optimum ****" << endl
		     << " ********************************************" << endl
		     << " solver returned: " << msg << endl;
	}

	//////////////////////
	// compute gradient //
	//////////////////////

	// phiCu = phi(x,y) = X^T*y
	_data->XTMultW(y, grad); // use grad in-place
	// phiGt = phi(x,y') = X^T*y'
	TheMatrix phiGt(_numFeatures, 1, SML::DENSE);
	_data->XTMultW(_data->labels(), phiGt);

	// gradient = phiCu - phiGt
	grad.Minus(phiGt);
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
