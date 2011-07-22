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
#include "timer.hpp"
#include "softmarginloss.hpp"
#include "configuration.hpp"
#include "costfactory.hpp"

using namespace std;

SoftMarginLoss::SoftMarginLoss(CModel* model, CConsVecData* data) :
	_data(data),
	_model(model),
	_numFeatures(_data->dim()),
	_numVariables(_data->NumOfLabel()),
	_numEqConstraints(_data->NumOfEqualities()),
	_numIneqConstraints(_data->NumOfInequalities()),
	_solver(_numVariables, _numEqConstraints, _numIneqConstraints),
	_costFunction(0),
	_costFactor(1.0),
	_gamma(1.0),
	_linearCostContribution(_numVariables, 1, SML::DENSE),
	_verbosity(0) {

	// get configurations
	Configuration &config = Configuration::GetInstance();

	if (config.IsSet("Loss.verbosity"))
		_verbosity = config.GetInt("Loss.verbosity");

	if (config.IsSet("Loss.gamma"))
		_gamma = config.GetDouble("Loss.gamma");

	if (_verbosity > 0) {
		cout << "[SoftMarginLoss] initialising..." << endl;

		cout << "[SoftMarginLoss] dataset consists of " << _numFeatures
			 << " features and " << _numVariables << " variables." << endl;
	}

	// get cost function
	_costFunction = CostFactory::GetCost(data);

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

SoftMarginLoss::~SoftMarginLoss() {

	if (_costFunction)
		delete _costFunction;
}

void
SoftMarginLoss::ComputeLoss(double& loss) {

	if (_verbosity > 1)
		cout << "[SoftMarginLoss::ComputeLoss] ..." << endl;

	TheMatrix grad(_numFeatures, 1, SML::DENSE);

	ComputeLossAndGradient(loss, grad);
}

void
SoftMarginLoss::ComputeLossAndGradient(double& loss, TheMatrix& grad) {

	CTimer totalTime;
	totalTime.Start();

	if (_verbosity > 1)
		cout << "[SoftMarginLoss::ComputeLossAndGradient] ..." << endl;


	/////////////////////////////////////////////////////////////////
	// compute coefficients f and constant part c of the objective //
	/////////////////////////////////////////////////////////////////

	CTimer computeCoefficientsTime;
	computeCoefficientsTime.Start();

	// the current weights
	TheMatrix& w = _model->GetW();

	if (_verbosity > 1) {

		cout << "[SoftMarginLoss::ComputeLossAndGradient] "
			 << "current weights: ";
		w.Print();
	}

	// f = Xw (intermediate result)
	TheMatrix f(_numVariables, 1);
	_data->XMultW(w, f);

	// c = c1 - c2 = const_cost - <y',Xw>

	// c1 = const_cost
	double c1 = _constantCostContribution;

	// c2 = <y',Xw>
	double c2;
	_data->labels().Dot(f, c2);

	double c = c1 - c2;

	// f = Xw + (1-2y')*gamma
	f.Add(_linearCostContribution);

	if (_verbosity > 2) {

		cout << "[SoftMarginLoss::ComputeLossAndGradient] "
			 << "constant term of objective: " << c << endl;
		cout << "[SoftMarginLoss::ComputeLossAndGradient] "
			 << "coefficients of objective: ";
		f.Print();
	}

	computeCoefficientsTime.Stop();

	CTimer setCoefficientsTime;
	setCoefficientsTime.Start();

	// set objective in linear solver
	_solver.SetObjective(f, c, LinearProgramSolver::MAXIMIZE);

	setCoefficientsTime.Stop();

	if (_verbosity > 1) {

		cout << "[SoftMarginLoss::ComputeLossAndGradient] "
		     << "time computing coefficients:"
		     << computeCoefficientsTime.WallclockTotal() << endl;
		cout << "[SoftMarginLoss::ComputeLossAndGradient] "
		     << "time setting coefficients  :"
		     << setCoefficientsTime.WallclockTotal() << endl;
	}

	///////////////////
	// solve the ILP //
	///////////////////

	CTimer solveTime;
	solveTime.Start();

	// the solution vector
	TheMatrix y(_numVariables, 1, SML::SPARSE);

	// the return message of the solver
	string msg;

	bool success = _solver.Solve(y, loss, msg);

	solveTime.Stop();

	if (_verbosity > 2) {

		cout << "[SoftMarginLoss::ComputeLossAndGradient] "
			 << "most offending solution: ";
		y.Print();
	}

	if (_verbosity > 1) {

		double costTerm;
		_linearCostContribution.Dot(y, costTerm);
		costTerm += _constantCostContribution;

		cout << "[SoftMarginLoss::ComputeLossAndGradient] "
		     << "loss Γ(y,y')<w,Θ(x,y')-Θ(x,y)> + Δ(y,y'): "
		     << (loss - costTerm)
		     << " + " << costTerm
		     << " = " << loss << endl;
	}

	if (!success) {

		cout << "[SoftMarginLoss::ComputeLossAndGradient]" << endl
		     << " **************** WARNING *******************" << endl
		     << " **** linear solver did not find optimum ****" << endl
		     << " ********************************************" << endl
		     << " solver returned: " << msg << endl;
	}

	if (_verbosity > 1)
		cout << "[SoftMarginLoss::ComputeLossAndGradient] "
		     << "time solving problem       :"
		     << solveTime.WallclockTotal() << endl;

	//////////////////////
	// compute gradient //
	//////////////////////

	CTimer computeGradientTime;
	computeGradientTime.Start();

	// phiCu = phi(x,y) = X^T*y
	_data->XTMultW(y, grad); // use grad in-place

	if (_verbosity > 2) {
		cout << "[SoftMarginLoss::ComputeLossAndGradient] "
			 << "phiCu: ";
		grad.Print();
	}

	// phiGt = phi(x,y') = X^T*y'
	TheMatrix phiGt(_numFeatures, 1, SML::DENSE);
	_data->XTMultW(_data->labels(), phiGt);

	if (_verbosity > 2) {

		cout << "[SoftMarginLoss::ComputeLossAndGradient] "
			 << "phiGt: ";
		phiGt.Print();
	}

	// gradient = phiCu - phiGt
	grad.Minus(phiGt);

	computeGradientTime.Stop();

	if (_verbosity > 1) {

		cout << "[SoftMarginLoss::ComputeLossAndGradient] "
			 << "gradient: ";
		grad.Print();

		cout << "[SoftMarginLoss::ComputeLossAndGradient] "
		     << "time computing gradient    :"
		     << computeGradientTime.WallclockTotal() << endl;
		cout << "[SoftMarginLoss::ComputeLossAndGradient] "
		     << "total time                 :"
		     << totalTime.WallclockTotal() << endl;
	}
}

void
SoftMarginLoss::ComputeCostContribution(const TheMatrix& groundTruth) {

	if (_verbosity > 2) {

		cout << "[SoftMarginLoss::ComputeCostContribution] "
			 << "computing loss contributions for " << endl;
		groundTruth.Print();
	}

	_costFunction->constantContribution(groundTruth, _constantCostContribution);
	_costFunction->linearContribution(groundTruth, _linearCostContribution);

	if (_verbosity > 2) {

		cout << "[SoftMarginLoss::ComputeCostContribution] "
			 << "constant contribution is " << _constantCostContribution
			 << endl;

		cout << "[SoftMarginLoss::ComputeCostContribution] "
			 << "linear contribution is " << endl;
		_linearCostContribution.Print();
	}
}
