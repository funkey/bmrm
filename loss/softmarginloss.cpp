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
	_numAuxiliaryVariables(0),
	_numAuxiliaryEqualities(0),
	_numAuxiliaryInequalities(0),
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

	// for linear gamma functions, we need to add auxiliary variables for the
	// resulting quadratic expression
	if (!_gammaConst)
		AddAuxiliaryVariables();

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

	//////////////////////////////////////////////
	// compute constant part c of the objective //
	//////////////////////////////////////////////

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

	// update constant margin contribution
	// m_c = -<Xw,y> = -<m_l,y>
	_m_l.Dot(_y, _m_c);
	_m_c *= -1;

	// constant part of objective
	double c = _g_c*_m_c + _c_c;

	//////////////////////////////////////////////////
	// compute linear coefficients of the objective //
	//////////////////////////////////////////////////

	// f is a vector of linear coefficients for the problem variables y' and the
	// auxiliary variables a:
	//
	// f = (f_y f_a)^T
	//
	// f_y = g_c*m_l + m_c*g_l + c_l + q_l
	// f_a = q_q

	TheMatrix f_y(_m_l);

	// _g_c*_m_l ...
	f_y.Scale(_g_c);

	// ... + m_c*g_l ... (only if g_l ≠ 0)
	if (!_gammaConst) {

		TheMatrix t(_g_l);
		t.Scale(_m_c);

		f_y.Add(t);
	}

	// ... + c_l
	f_y.Add(_c_l);

	if (_verbosity > 2) {

		cout << "[SoftMarginLoss::ComputeLossAndGradient] "
			 << "constant term of objective: " << c << endl;
		cout << "[SoftMarginLoss::ComputeLossAndGradient] "
			 << "linear coefficients of objective: ";
		f_y.Print();
	}

	////////////////////////////////////////////////////////////////////////
	// compute auxiliary coefficients for quadratic term of the objective //
	////////////////////////////////////////////////////////////////////////

	CTimer setCoefficientsTime;

	if (!_gammaConst) {

		TheMatrix f_a(_numAuxiliaryVariables, 1);

		int varnum = 0;

		// f_a is vectorized upper triangle of Q ≈ Γ_l*M_l^T, where Q is lower
		// triabgle zero matrix such that x^T*Q*x = x^T*(Γ_l*M_l^T)*x for every
		// x.
		//
		// f_a = vec(_g_l*_m_l^T)
		//
		for (int i = 0; i < _numVariables - _numAuxiliaryVariables; i++) {

			double g_l_i;
			double m_l_i;

			_g_l.Get(i, g_l_i);
			_m_l.Get(i, m_l_i);

			for (int j = i + 1; j < _numVariables - _numAuxiliaryVariables; j++) {

				double g_l_j;
				double m_l_j;

				_g_l.Get(j, g_l_j);
				_m_l.Get(j, m_l_j);

				f_a.Set(varnum, g_l_i*m_l_j + g_l_j*m_l_i);
				varnum += 3;
			}

			// diagonal of Q is linear for binary y, therefore we add it to f_y
			double current;
			f_y.Get(i, current);
			f_y.Set(i, current + g_l_i*m_l_i);
		}

		if (_verbosity > 2) {

			cout << "[SoftMarginLoss::ComputeLossAndGradient] "
				 << "augmented linear coefficients of objective: ";
			f_y.Print();

			cout << "[SoftMarginLoss::ComputeLossAndGradient] "
				 << "linearized quadratic coefficients of objective: ";
			f_a.Print();
		}

		// join linear and linearized quadratic coefficients
		// TODO: reuse f
		TheMatrix f(_numVariables, 1);

		for (int i = 0; i < _numVariables - _numAuxiliaryVariables; i++) {
			double value;
			f_y.Get(i, value);
			f.Set(i, value);
		}
		for (int i = _numVariables - _numAuxiliaryVariables; i < _numVariables; i++) {
			double value;
			f_a.Get(i - (_numVariables - _numAuxiliaryVariables), value);
			f.Set(i, value);
		}

		computeCoefficientsTime.Stop();
		setCoefficientsTime.Start();

		// set objective in linear solver
		_solver->SetObjective(f, c, LinearProgramSolver::MAXIMIZE);

		setCoefficientsTime.Stop();

	} else {

		computeCoefficientsTime.Stop();
		setCoefficientsTime.Start();

		// set objective in linear solver
		_solver->SetObjective(f_y, c, LinearProgramSolver::MAXIMIZE);

		setCoefficientsTime.Stop();
	}


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
	// TODO: reuse y_
	TheMatrix y_(_numVariables - _numAuxiliaryVariables, 1, SML::SPARSE);

	// the return message of the solver
	string msg;
	bool success;

	if (!_gammaConst) {

		TheMatrix y_all(_numVariables, 1, SML::SPARSE);
		success = _solver->Solve(y_all, loss, msg);

		// read back the problem variable part only
		for (int i = 0; i < _numVariables - _numAuxiliaryVariables; i++) {
			double value;
			y_all.Get(i, value);
			y_.Set(i, value);
		}

	} else {

		success = _solver->Solve(y_, loss, msg);
	}

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

	// gamma = g_c + <g_l,y'>
	double gamma;
	_g_l.Dot(y_, gamma);
	gamma += _g_c;

	// gradient = gamma*(phiCu - phiGt)
	grad.Minus(phiGt);
	grad.Scale(gamma);

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
					"unkown Loss.gammaFunctionType value",
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
SoftMarginLoss::AddAuxiliaryVariables() {

	/*
	 * Assuming we have a quadratic part in our objective of the form y^T*Q*y,
	 * we introduce one auxiliary variable for each pair (i,j), i ≠ j, i < j:
	 *
	 * y^T*Q*y = y^T*Q'*y,
	 *
	 * where Q' is an upper triangular matrix with Q'_ij = Q_ij + Q_ji (i < j)
	 * and Q'_ii = Q_ii.
	 *
	 * The diagonal entries of Q' contribute linearly to the objective (since y
	 * is binary). It remains to add auxiliary variables for the upper triangle
	 * without the diagonal, that is, for n(n-1)/2 entries.
	 */

	if (_verbosity > 1)
		cout << "[SoftMarginLoss::AddAuxiliaryVariables] "
		     << "computing number of auxiliary variables for "
			 << _numVariables << " problem variables" << endl;

	int numEntries =
			_numVariables*(_numVariables-1)/2;

	if (_verbosity > 2)
		cout << "[SoftMarginLoss::AddAuxiliaryVariables] Q has "
		     << numEntries << " entries in upper triangle" << endl;

	// a_ij, a_-ij, a_i-j
	_numAuxiliaryVariables = numEntries*3;

	_numAuxiliaryInequalities = numEntries;
	_numAuxiliaryEqualities   = numEntries*2;

	if (_verbosity > 2)
		cout << "[SoftMarginLoss::AddAuxiliaryVariables] add "
		     << _numAuxiliaryVariables << " auxiliary variables, "
		     << _numAuxiliaryEqualities << " auxiliary equalities, "
		     << _numAuxiliaryInequalities << " auxiliary inequalities"
		     << endl;

	_numVariables    += _numAuxiliaryVariables;
	_numEqualities   += _numAuxiliaryEqualities;
	_numInequalities += _numAuxiliaryInequalities;
}

void
SoftMarginLoss::SetLinearConstraints() {

	if (_numEqualities > 0) {

		// check if ground-truth fulfills the constraints
		if (_verbosity > 0)
			if (_numEqualities - _numAuxiliaryEqualities > 0)
				CheckSolutionIntegrity(
						_data->GetEqualityCoefs(),
						_y, 0,
						_data->GetEqualityValues());

		if (_gammaConst) { // use equalities from data directly

			_solver->SetEqualities(
					_data->GetEqualityCoefs(),
					_data->GetEqualityValues());

		} else {

			if (_verbosity > 1)
				cout << "[SoftMarginLoss::SetLinearConstraints] "
				     << "creating equalities coefficients and values for "
				     << _numEqualities << " equalities and "
				     << _numVariables << " variables"
				     << endl;

			// create equality coefficent matrix and value vector
			TheMatrix eqCoefs(_numEqualities, _numVariables, SML::SPARSE);
			TheMatrix eqValues(_numEqualities, 1, SML::DENSE);

			if (_verbosity > 1)
				cout << "[SoftMarginLoss::SetLinearConstraints] "
				     << "filling coefficients and values..."
				     << endl;

			// equalities from the data
			for (int j = 0; j < _numEqualities - _numAuxiliaryEqualities; j++) {
				for (int i = 0; i < _numVariables - _numAuxiliaryVariables; i++) {

					double coef;
					_data->GetEqualityCoefs().Get(j, i, coef);

					eqCoefs.Set(j, i, coef);
				}

				double value;
				_data->GetEqualityValues().Get(j, value);
				eqValues.Set(j, value);
			}

			if (_verbosity > 1)
				cout << "[SoftMarginLoss::SetLinearConstraints] "
				     << "...done with data part..."
				     << endl;

			// equalities from the auxiliary variables
			int auxEqNum  = _numEqualities - _numAuxiliaryEqualities;
			int auxVarNum = _numVariables  - _numAuxiliaryVariables;

			for (int i = 0; i < _numVariables - _numAuxiliaryVariables; i++) {
				for (int j = i+1; j < _numVariables - _numAuxiliaryVariables; j++) {

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

			if (_verbosity > 1)
				cout << "[SoftMarginLoss::SetLinearConstraints] "
				     << "...done with auxiliary part."
				     << endl;

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

		// check if ground-truth fulfills the constraints
		if (_verbosity > 0)
			if (_numInequalities - _numAuxiliaryInequalities > 0)
				CheckSolutionIntegrity(
						_data->GetInequalityCoefs(),
						_y, -1,
						_data->GetInequalityValues());

		if (_gammaConst) { // use inequalities from data directly

			_solver->SetInequalities(
					_data->GetInequalityCoefs(),
					_data->GetInequalityValues());

		} else {

			if (_verbosity > 1)
				cout << "[SoftMarginLoss::SetLinearConstraints] "
				     << "creating inequalities coefficients and values for "
				     << _numInequalities << " inequalities and "
				     << _numVariables << " variables"
				     << endl;

			// create equality coefficent matrix and value vector
			TheMatrix ineqCoefs(_numInequalities, _numVariables, SML::SPARSE);
			TheMatrix ineqValues(_numInequalities, 1, SML::DENSE);

			if (_verbosity > 1)
				cout << "[SoftMarginLoss::SetLinearConstraints] "
				     << "filling coefficients and values..."
				     << endl;

			// inequalities from the data
			for (int j = 0; j < _numInequalities - _numAuxiliaryInequalities; j++) {
				for (int i = 0; i < _numVariables - _numAuxiliaryVariables; i++) {

					double coef;
					_data->GetInequalityCoefs().Get(j, i, coef);

					ineqCoefs.Set(j, i, coef);
				}

				double value;
				_data->GetInequalityValues().Get(j, value);
				ineqValues.Set(j, value);
			}

			if (_verbosity > 1)
				cout << "[SoftMarginLoss::SetLinearConstraints] "
				     << "...done with data part..."
				     << endl;

			// inequalities from the auxiliary variables
			int auxIneqNum  = _numInequalities - _numAuxiliaryInequalities;
			int auxVarNum   = _numVariables  - _numAuxiliaryVariables;

			for (int i = 0; i < _numVariables - _numAuxiliaryVariables; i++) {
				for (int j = i+1; j < _numVariables - _numAuxiliaryVariables; j++) {

					// 1*a_ij + 1*a_i-j + 1*a_-ij <= 1
					ineqCoefs.Set(auxIneqNum, auxVarNum,     1.0); // a_ij
					ineqCoefs.Set(auxIneqNum, auxVarNum + 1, 1.0); // a_i-j
					ineqCoefs.Set(auxIneqNum, auxVarNum + 2, 1.0); // a_-ij
					ineqValues.Set(auxIneqNum, 1);

					auxIneqNum++;
					auxVarNum += 3;
				}
			}

			if (_verbosity > 1)
				cout << "[SoftMarginLoss::SetLinearConstraints] "
				     << "...done with auxiliary part."
				     << endl;

			if (auxIneqNum != _numInequalities)
				cout << "[SoftMarginLoss::SetLinearConstraints] " << auxIneqNum
				     << " linear inequalities set instead of " << _numInequalities
				     << endl;

			if (_verbosity > 2) {

				cout << "[SoftMarginLoss::SetLinearConstraints] "
				     << "linear inequality constraints:" << endl;
				cout << endl;
				ineqCoefs.Print();
				cout << endl;
				ineqValues.Print();
			}

			// set them
			_solver->SetInequalities(ineqCoefs, ineqValues);
		}
	}
}

void
SoftMarginLoss::CheckSolutionIntegrity(
		const TheMatrix& A,
		const TheMatrix& y,
		int relation,
		const TheMatrix& b) {

	TheMatrix x(b.Rows(), 1);
	A.Dot(y, x);

	bool failed = false;

	for (int i = 0; i < A.Rows(); i++) {

		double x_i;
		double b_i;

		x.Get(i, x_i);
		b.Get(i, b_i);

		switch (relation) {

			case -1:
				if (x_i > b_i) {
					cout << "[SoftMarginLoss::CheckSolutionIntegrity] "
					     << "ground-truth violates " << i << "th "
					     << "\"<=\" constraint! (" << x_i << " > "
					     << b_i << ")" << endl;
					failed = true;
				}
				break;

			case 0:
				if (x_i != b_i) {
					cout << "[SoftMarginLoss::CheckSolutionIntegrity] "
					     << "ground-truth violates " << i << "th "
					     << "\"==\" constraint! (" << x_i << " != "
					     << b_i << ")" << endl;
					failed = true;
				}
				break;

			case 1:
				if (x_i < b_i) {
					cout << "[SoftMarginLoss::CheckSolutionIntegrity] "
					     << "ground-truth violates " << i << "th "
					     << "\">=\" constraint! (" << x_i << " < "
					     << b_i << ")" << endl;
					failed = true;
				}
				break;
		}
	}

	if (!failed)
		cout << "[SoftMarginLoss::CheckSolutionIntegrity] "
			 << "ground-truth satisfies "
			 << (relation == 0 ? "" : "in") << "equality constraints"
			 << endl;
}
