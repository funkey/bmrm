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
 * Created: (07/18/2011)
 */

#include <sstream>

#include <CplexSolver.hpp>

#include "timer.hpp"
#include "configuration.hpp"

CplexSolver::CplexSolver(
		unsigned int numVariables,
		unsigned int numEqConstraints,
		unsigned int numIneqConstraints) :
	LinearProgramSolver(
			numVariables,
			numEqConstraints,
			numIneqConstraints),
	_variables(_env, numVariables, 0, 1, ILOINT),
	_objective(_env),
	_coefs(_env, numVariables),
	_model(_env),
	_cplex(_model),
	_verbosity(0) {

	// get configurations
	Configuration &config = Configuration::GetInstance();

	if (config.IsSet("Cplex.verbosity"))
		_verbosity = config.GetInt("Cplex.verbosity");

	// setup Ilo environment
	if (_verbosity == 0)
		_env.setOut(_env.getNullStream());

		// add all variables to the model
	for (int i = 0; i < numVariables; i++)
		_model.add(
				IloRange(
						_env,
						_variables[i].getLB(),
						_variables[i],
						_variables[i].getUB()));

	if (_verbosity > 1)
		cout << "[CplexSolver] created " << numVariables
			 << " binary variables" << endl;

	// add the objective to the model (it will be filled later)
	_model.add(_objective);
}

CplexSolver::~CplexSolver() {

	// destroy Ilo environment
	_env.end();
}

void
CplexSolver::SetObjective(const TheMatrix& a, double c, Sense sense) {

	// TODO: check size of a

	// normalize the values of a and c to an intervall between -1 and 1
	double max = c;
	for (unsigned int i = 0; i < _numVariables; i++) {

		double value;
		a.Get(i, value);

		if (abs(value) > max)
			max = abs(value);
	}
	_scale = 1.0/max;

	// set sense of objective
	if (sense == MINIMIZE)
		_objective.setSense(IloObjective::Minimize);
	else
		_objective.setSense(IloObjective::Maximize);

	// set the constant value of the objective
	_objective.setConstant(static_cast<IloNum>(_scale*c));

	// set the coefficients for all variables
	for (unsigned int i = 0; i < _numVariables; i++) {

		double value;
		a.Get(i, value);

		_coefs[i] = static_cast<IloNum>(_scale*value);
	}

	_objective.setLinearCoefs(_variables, _coefs);

	if (_verbosity > 2)
		cout << "[CplexSolver] model after setting objective: "
			 << endl << _model << endl;
}

void
CplexSolver::SetEqualities(const TheMatrix& A, const TheMatrix& b) {

	if (_verbosity > 1)
		cout << "[CplexSolver::SetEqualities] setting "
		     << _numEqConstraints
		     << " equality constraints" << endl;

	CTimer totalTimer;
	totalTimer.Start();

	// 0 --> equality constraints
	SetConstraints(A, b, 0);

	totalTimer.Stop();

	if (_verbosity > 1)
		cout << "[CplexSolver::SetEqualities] "
		     << _numEqConstraints << " set in "
		     << totalTimer.WallclockTotal()
		     << " seconds" << endl;

	if (_verbosity > 2)
		cout << "[CplexSolver::SetEqualities] "
		     << "model after setting equality constraints: "
			 << endl << _model << endl;
}

void
CplexSolver::SetInequalities(const TheMatrix& A, const TheMatrix& b) {

	if (_verbosity > 1)
		cout << "[CplexSolver::SetInequalities] setting "
		     << _numIneqConstraints
		     << " inequality constraints" << endl;

	CTimer totalTimer;
	totalTimer.Start();

	// -1 -->  "<=" inequality constraints
	SetConstraints(A, b, -1);

	totalTimer.Stop();

	if (_verbosity > 1)
		cout << "[CplexSolver::SetInequalities] "
		     << _numIneqConstraints << " set in "
		     << totalTimer.WallclockTotal()
		     << " seconds" << endl;

	if (_verbosity > 2)
		cout << "[CplexSolver::SetInequalities] "
		     << "model after setting inequality constraints: "
			 << endl << _model << endl;
}

void CplexSolver::SetConstraints(const TheMatrix& A, const TheMatrix& b, int relation) {

	IloRangeArray constraints(_env);
	unsigned int  numConstraints = A.Rows();

	CTimer initTimer;
	initTimer.Start();

	unsigned int  numEntries = 0;
	double*       coefs      = new double[_numVariables];
	unsigned int* indices    = new unsigned int[_numVariables];

	initTimer.Stop();

	if (_verbosity > 1)
		cout << "[CplexSolver::SetConstraints] "
		     << "time initialising           : "
		     << initTimer.WallclockTotal() << endl;

	CTimer extractTimer;
	CTimer setTimer;

	for (unsigned int j = 0; j < numConstraints; j++) {

		if (_verbosity > 1 && j > 0) {
			if (j % 1000 == 0) {

				cout << "[CplexSolver::SetConstraints] "
				     << j << " constraints set so far"
				     << " (mean extraction time: "
				     << extractTimer.WallclockTotal()/j
				     << ", mean set time: "
				     << setTimer.WallclockTotal()/j
				     << ")" << endl;
			}
		}

		extractTimer.Start();

		// get the whole row
		A.GetRow(j, numEntries, coefs, indices);

		// get the rhs value
		double value;
		b.Get(j, value);

		extractTimer.Stop();

		setTimer.Start();

		// set the bounds
		IloRange constraint(
				_env,
				(relation < 0 ? -IloInfinity : value),
				(relation > 0 ?  IloInfinity : value));

		// set the coefficients
		for (unsigned int i = 0; i < numEntries; i++) {

			constraint.setLinearCoef(
					_variables[indices[i]],
					static_cast<IloNum>(coefs[i]));
		}

		// add to the constraint array
		constraints.add(constraint);

		setTimer.Stop();
	}

	delete[] coefs;
	delete[] indices;

	_model.add(constraints);
}

bool
CplexSolver::Solve(TheMatrix& x, double& value, string& msg) {

	_cplex.solve();

	// get solver result message
	stringstream ss;
	ss << _cplex.getStatus() << flush;
	msg = ss.str();

	if (_cplex.getStatus() != IloAlgorithm::Optimal)
		return false;

	// extract solution
	IloNumArray values(_env);
	_cplex.getValues(values, _variables);
	for (unsigned int i = 0; i < _numVariables; i++)
		x.Set(i, static_cast<double>(values[i]));

	// get current value of the objective
	value = static_cast<double>(_cplex.getObjValue())/_scale;

	return true;
}
