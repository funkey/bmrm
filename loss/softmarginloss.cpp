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
	_numEqualities(_data->NumOfEqualities()),
	_numInequalities(_data->NumOfInequalities()),
	_numAuxilaryVariables(0),
	_numAuxilaryEqualities(0),
	_numAuxilaryInequalities(0),
	_solver(0),
	_costFunction(0),
	_costFactor(1.0),
	_c_l(_numVariables, 1, SML::DENSE),
	_g_l(_numVariables, 1, SML::DENSE),
	_m_l(_numVariables, 1, SML::DENSE),
	_y(_data->labels()),
	_verbosity(0) {

	// get configurations
	Configuration &config = Configuration::GetInstance();

	if (config.IsSet("Loss.verbosity"))
		_verbosity = config.GetInt("Loss.verbosity");

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

	// calculate cost contribution (does not change anymore)
	ComputeCostContribution(_y);

	// calculate gamma contribution (does not change anymore)
	ComputeGammaContribution(_y);

	// for linear gamma functions, we need to add auxilary variables for the
	// resulting quadratic expression
	if (!_gammaConst)
		AddAuxilaryVariables();

	_solver = new CplexSolver(
			_numVariables,
			_numEqualities,
			_numInequalities);

	// set linear constraints (do not change anymore)
	SetLinearConstraints();
}

SoftMarginLoss::~SoftMarginLoss() {

	if (_costFunction)
		delete _costFunction;

	if (_solver)
		delete _solver;
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

	// the current weights
	TheMatrix& w = _model->GetW();

	////////////////////////////////////////////////////////////////////////
	// compute linear coefficients f and constant part c of the objective //
	////////////////////////////////////////////////////////////////////////

	CTimer computeCoefficientsTime;
	computeCoefficientsTime.Start();

	if (_verbosity > 1) {

		cout << "[SoftMarginLoss::ComputeLossAndGradient] "
			 << "current weights: ";
		w.Print();
	}

	// update linear margin contribution
	// m_l = Xw
	_data->XMultW(w, _m_l);

	// c = g_c*m_c + c_c

	// update constant margin contribution
	// m_c = -<Xw,y> = -<m_l,y>
	_m_l.Dot(_y, _m_c);
	_m_c *= -1;

	double c = _g_c*_m_c + _c_c;

	// f = g_c*m_l + m_c*g_l + c_l

	// _g_c*_m_l ...
	TheMatrix f(_m_l);
	f.Scale(_g_c);

	// ... + m_c*g_l ...
	if (!_gammaConst) {

		TheMatrix t(_g_l);
		t.Scale(_m_c);
		f.Add(t);
	}

	// ... + c_l
	f.Add(_c_l);

	if (_verbosity > 2) {

		cout << "[SoftMarginLoss::ComputeLossAndGradient] "
			 << "constant term of objective: " << c << endl;
		cout << "[SoftMarginLoss::ComputeLossAndGradient] "
			 << "coefficients of objective: ";
		f.Print();
	}

	/////////////////////////////////////////////
	// compute quadratic term of the objective //
	/////////////////////////////////////////////

	if (!_gammaConst) {

		TheMatrix q(_numVariables, _numVariables);

		for (int i = 0; i < _numVariables; i++)
			for (int j = 0; j < _numVariables; j++) {

				double g_l_i;
				double m_l_j;

				_g_l.Get(i, g_l_i);
				_m_l.Get(j, m_l_j);

				q.Set(i, j, g_l_i*m_l_j);
			}

		if (_verbosity > 2) {

			cout << "[SoftMarginLoss::ComputeLossAndGradient] "
				 << "quadradic term of objective: " << endl;
			q.Print();
		}
	}

	computeCoefficientsTime.Stop();

	CTimer setCoefficientsTime;
	setCoefficientsTime.Start();

	// set objective in linear solver
	_solver->SetObjective(f, c, LinearProgramSolver::MAXIMIZE);

	// TODO:
	// set quadratic term

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
	TheMatrix y_(_numVariables, 1, SML::SPARSE);

	// the return message of the solver
	string msg;

	bool success = _solver->Solve(y_, loss, msg);

	solveTime.Stop();

	if (_verbosity > 2) {

		cout << "[SoftMarginLoss::ComputeLossAndGradient] "
			 << "most offending solution: ";
		y_.Print();
	}

	if (_verbosity > 1) {

		double costTerm;
		_c_l.Dot(y_, costTerm);
		costTerm += _c_c;

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

	// phiCu = phi(x,y') = X^T*y'
	_data->XTMultW(y_, grad); // use grad in-place

	if (_verbosity > 2) {
		cout << "[SoftMarginLoss::ComputeLossAndGradient] "
			 << "phiCu: ";
		grad.Print();
	}

	// phiGt = phi(x,y) = X^T*y
	TheMatrix phiGt(_numFeatures, 1, SML::DENSE);
	_data->XTMultW(_y, phiGt);

	if (_verbosity > 2) {

		cout << "[SoftMarginLoss::ComputeLossAndGradient] "
			 << "phiGt: ";
		phiGt.Print();
	}

	// gradient = gamma*(phiCu - phiGt)
	grad.Minus(phiGt);
	grad.Scale(_g_c);

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

	_costFunction->constantContribution(groundTruth, _c_c);
	_costFunction->linearContribution(groundTruth, _c_l);

	if (_verbosity > 2) {

		cout << "[SoftMarginLoss::ComputeCostContribution] "
			 << "constant contribution is " << _c_c
			 << endl;

		cout << "[SoftMarginLoss::ComputeCostContribution] "
			 << "linear contribution is " << endl;
		_c_l.Print();
	}
}

void
SoftMarginLoss::ComputeGammaContribution(const TheMatrix& groundTruth) {

	if (_verbosity > 2) {

		cout << "[SoftMarginLoss::ComputeGammaContribution] "
			 << "computing gamma contributions for " << endl;
		groundTruth.Print();
	}

	// default: constant gamma
	_g_c        = 1.0;
	_gammaConst = true;

	Configuration &config = Configuration::GetInstance();

	if (config.IsSet("Loss.gammaFunctionType")) {

		string gammaFunctionType = config.GetString("Loss.gammaFunctionType");

		if (gammaFunctionType == "CONSTANT") {

			if (config.IsSet("Loss.gammaConstant"))
				_g_c = config.GetDouble("Loss.gammaConstant");

		} else if (gammaFunctionType == "COST_FUNCTION") {

			_g_c = _c_c;
			_g_l = _c_l;

			_gammaConst = false;

		} else
			throw CBMRMException(
					"unkown Loss.gammaConstant value",
					"SoftMarginLoss::ComputeGammaContribution");
	}

	if (_verbosity > 2) {

		cout << "[SoftMarginLoss::ComputeGammaContribution] "
			 << "constant contribution is " << _g_c
			 << endl;

		cout << "[SoftMarginLoss::ComputeGammaContribution] "
			 << "linear contribution is " << endl;

		_g_l.Print();
	}
}

void
SoftMarginLoss::AddAuxilaryVariables() {

	/*
	 * Assuming we have a quadratic part in our objective of the form y^T*Q*y,
	 * we introduce one auxilary variable for each pair (i,j), i ≠ j, i < j:
	 *
	 * y^T*Q*y = y^T*Q'*y,
	 *
	 * where Q' is an upper triangular matrix with Q'_ij = Q_ij + Q_ji (i < j)
	 * and Q'_ii = Q_ii.
	 *
	 * The diagonal entries of Q' contribute linearly to the objective (since y
	 * is binary). It remains to add auxilary variables for the upper triangle
	 * without the diagonal, that is, for n(n-1)/2 entries.
	 */

	int numEntries =
			_numVariables*(_numVariables-1)/2;

	// a_ij, a_-ij, a_i-j
	_numAuxilaryVariables = numEntries*3;

	_numAuxilaryInequalities = numEntries;
	_numAuxilaryEqualities   = numEntries*2;

	_numVariables    += _numAuxilaryVariables;
	_numEqualities   += _numAuxilaryEqualities;
	_numInequalities += _numAuxilaryInequalities;
}

void
SoftMarginLoss::SetLinearConstraints() {

	if (_numEqualities > 0) {

		if (_gammaConst) // use equalities from data directly

			_solver->SetEqualities(
					_data->GetEqualityCoefs(),
					_data->GetEqualityValues());
		else {

			// create equality coefficent matrix and value vector
			TheMatrix eqCoefs(_numEqualities, _numVariables, SML::SPARSE);
			TheMatrix eqValues(_numEqualities, SML::DENSE);

			// equalities from the data
			for (int i = 0; i < _numVariables - _numAuxilaryVariables; i++) {
				for (int j = 0; j < _numEqualities - _numAuxilaryEqualities; j++) {

					double coef;
					_data->GetEqualityCoefs().Get(j, i, coef);

					eqCoefs.Set(j, i, coef);
				}
			}

			// equalities from the auxilary variables
			int auxEqNum  = _numEqualities - _numAuxilaryEqualities;
			int auxVarNum = _numVariables;

			for (int i = 0; i < _numVariables; i++) {
				for (int j = i+1; j < _numVariables; j++) {

					// 1*a_ij + 1*a_i-j - 1*y_i == 0
					eqCoefs.Set(auxEqNum, auxVarNum,     1.0); // a_ij
					eqCoefs.Set(auxEqNum, auxVarNum + 1, 1.0); // a_i-j
					eqCoefs.Set(auxEqNum, i,            -1.0); // y_i
					eqValues.Set(auxEqNum, 0);

					auxEqNum++;

					// 1*a_ij + 1*a_-ij - 1*y_j == 0
					eqCoefs.Set(auxEqNum, auxVarNum,     1.0); // a_ij
					eqCoefs.Set(auxEqNum, auxVarNum + 2, 1.0); // a_-ij
					eqCoefs.Set(auxEqNum, j,            -1.0); // y_j
					eqValues.Set(auxEqNum, 0);

					auxEqNum++;

					auxVarNum += 3;
				}
			}

			if (auxEqNum != _numEqualities)
				cout << "[SoftMarginLoss::SetLinearConstraints] " << auxEqNum
				     << " linear equalities set instead of " << _numEqualities
				     << endl;

			if (_verbosity > 2) {

				cout << "[SoftMarginLoss::SetLinearConstraints] "
				     << "linear equality constraints:" << endl;
				cout << endl;
				eqCoefs.Print();
				cout << endl;
				eqValues.Print();
			}

			// set them
			_solver->SetEqualities(eqCoefs, eqValues);
		}
	}

	if (_numInequalities > 0) {

		if (_gammaConst) // use inequalities from data directly

			_solver->SetInequalities(
					_data->GetInequalityCoefs(),
					_data->GetInequalityValues());
		else {

			// create inequality coefficent matrix
			// create inequality value vector
			// set them
		}
	}
}
