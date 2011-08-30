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
CplexSolver::SetObjective(const TheMatrix& a, double constant, Sense sense) {

	// TODO: check size of a

	// set sense of objective
	if (sense == MINIMIZE)
		_objective.setSense(IloObjective::Minimize);
	else
		_objective.setSense(IloObjective::Maximize);

	// set the constant value of the objective
	_objective.setConstant(static_cast<IloNum>(constant));

	// set the coefficients for all variables
	for (unsigned int i = 0; i < _numVariables; i++) {

		double value;
		a.Get(i, value);

		_coefs[i] = static_cast<IloNum>(value);
	}

	_objective.setLinearCoefs(_variables, _coefs);

	if (_verbosity > 2)
		cout << "[CplexSolver] model after setting objective: "
			 << endl << _model << endl;
}

void
CplexSolver::SetEqualities(const TheMatrix& A, const TheMatrix& b) {

	IloRangeArray constraints(_env);

	if (_verbosity > 1)
		cout << "[CplexSolver::SetEqualities] setting "
		     << _numIneqConstraints
		     << " equality constraints" << endl;

	CTimer initTimer;
	initTimer.Start();

	unsigned int  numEntries = 0;
	double*       coefs      = new double[_numVariables];
	unsigned int* indices    = new unsigned int[_numVariables];

	initTimer.Stop();

	if (_verbosity > 1)
		cout << "[CplexSolver::SetEqualities] "
		     << "time initialising           : "
		     << initTimer.WallclockTotal() << endl;

	for (unsigned int j = 0; j < _numEqConstraints; j++) {

		if (_verbosity > 1)
			if (j % 100 == 0)
				cout << "[CplexSolver::SetEqualities] "
				     << j << " constraints set so far"
				     << endl;

		CTimer extractTimer;
		extractTimer.Start();

		// get the whole row
		A.GetRow(j, numEntries, coefs, indices);

		// get the rhs value
		double value;
		b.Get(j, value);

		extractTimer.Stop();

		CTimer setTimer;
		setTimer.Start();

		// set the bounds
		IloRange constraint(_env, value, value);

		// set the coefficients
		for (unsigned int i = 0; i < numEntries; i++) {

			constraint.setLinearCoef(
					_variables[indices[i]],
					static_cast<IloNum>(coefs[i]));
		}

		// add to the constraint array
		constraints.add(constraint);

		setTimer.Stop();

		if (_verbosity > 1) {

			cout << "[CplexSolver::SetEqualities] "
			     << "time extracting coefficients: "
			     << extractTimer.WallclockTotal() << endl;
			cout << "[CplexSolver::SetEqualities] "
			     << "time setting coefficients   : "
			     << setTimer.WallclockTotal() << endl;
		}
	}

	_model.add(constraints);

	if (_verbosity > 2)
		cout << "[CplexSolver::SetEqualities] "
		     << "model after setting equality constraints: "
			 << endl << _model << endl;
}

void
CplexSolver::SetInequalities(const TheMatrix& C, const TheMatrix& d) {

	IloRangeArray constraints(_env);

	if (_verbosity > 1)
		cout << "[CplexSolver::SetInequalities] setting "
		     << _numIneqConstraints
		     << " inequality constraints" << endl;

	CTimer initTimer;
	initTimer.Start();

	unsigned int  numEntries = 0;
	double*       coefs      = new double[_numVariables];
	unsigned int* indices    = new unsigned int[_numVariables];

	initTimer.Stop();

	if (_verbosity > 1)
		cout << "[CplexSolver::SetInequalities] "
		     << "time initialising           : "
		     << initTimer.WallclockTotal() << endl;

	for (unsigned int j = 0; j < _numIneqConstraints; j++) {

		if (_verbosity > 1)
			if (j % 100 == 0)
				cout << "[CplexSolver::SetInequalities] "
				     << j << " constraints set so far"
				     << endl;

		CTimer extractTimer;
		extractTimer.Start();

		// get the whole row
		C.GetRow(j, numEntries, coefs, indices);

		// get the rhs value
		double value;
		d.Get(j, value);

		extractTimer.Stop();

		CTimer setTimer;
		setTimer.Start();

		// set the bounds
		IloRange constraint(_env, -IloInfinity, value);

		// set the coefficients
		for (unsigned int i = 0; i < numEntries; i++) {

			constraint.setLinearCoef(
					_variables[indices[i]],
					static_cast<IloNum>(coefs[i]));
		}

		// add to the constraint array
		constraints.add(constraint);

		setTimer.Stop();

		if (_verbosity > 1) {

			cout << "[CplexSolver::SetInequalities] "
			     << "time extracting coefficients: "
			     << extractTimer.WallclockTotal() << endl;
			cout << "[CplexSolver::SetInequalities] "
			     << "time setting coefficients   : "
			     << setTimer.WallclockTotal() << endl;
		}
	}

	delete[] coefs;
	delete[] indices;

	_model.add(constraints);

	if (_verbosity > 2)
		cout << "[CplexSolver::SetInequalities] "
		     << "model after setting inequality constraints: "
			 << endl << _model << endl;
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
	value = static_cast<double>(_cplex.getObjValue());

	return true;
}
