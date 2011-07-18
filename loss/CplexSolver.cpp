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
	_model(_env),
	_cplex(_model) {

	// setup Ilo environment
	//_env.setOut(env.getNullStream());

		// add all variables to the model
	for (int i = 0; i < numVariables; i++)
		_model.add(
				IloRange(
						_env,
						_variables[i].getLB(),
						_variables[i],
						_variables[i].getUB()));
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
CplexSolver::SetObjective(const TheMatrix& a, double constant) {

	// TODO: check size of a

	// set the constant value of the objective
	_objective.setConstant(static_cast<IloNum>(constant));

	// set the coefficients for all variables
	for (unsigned int i = 0; i < _numVariables; i++) {

		double value;
		a.Get(i, value);

		_objective.setLinearCoef(_variables[i], static_cast<IloNum>(value));
	}

	cout << "[CplexSolver] model after setting objective: "
	     << endl << _model << endl;
}

void
CplexSolver::SetEqualities(const TheMatrix& A, const TheMatrix& b) {

	IloRangeArray constraints(_env);

	for (unsigned int j = 0; j < _numEqConstraints; j++) {

		// get the rhs value
		double value;
		b.Get(j, value);

		// set the bounds
		IloRange constraint(_env, value, value);

		// set the coefficients
		for (unsigned int i = 0; i < _numVariables; i++) {

			double coefficient;
			A.Get(j, i, coefficient);

			if (coefficient == 0.0)
				continue;

			constraint.setLinearCoef(
					_variables[i],
					static_cast<IloNum>(coefficient));
		}

		// add to the constraint array
		constraints.add(constraint);
	}

	_model.add(constraints);

	cout << "[CplexSolver] model after setting equality constraints: "
	     << endl << _model << endl;
}

void
CplexSolver::SetInequalities(const TheMatrix& C, const TheMatrix& d) {

	IloRangeArray constraints(_env);

	for (unsigned int j = 0; j < _numIneqConstraints; j++) {

		// get the rhs value
		double value;
		d.Get(j, value);

		// set the bounds
		IloRange constraint(_env, -IloInfinity, value);

		// set the coefficients
		for (unsigned int i = 0; i < _numVariables; i++) {

			double coefficient;
			C.Get(j, i, coefficient);

			if (coefficient == 0.0)
				continue;

			constraint.setLinearCoef(
					_variables[i],
					static_cast<IloNum>(coefficient));
		}

		// add to the constraint array
		constraints.add(constraint);
	}

	_model.add(constraints);

	cout << "[CplexSolver] model after setting inequality constraints: "
	     << endl << _model << endl;
}

bool
CplexSolver::Solve(TheMatrix& x, double& value, string& msg) {

	_cplex.solve();

	// get solver result message
	msg = "";
	stringstream ss(msg);
	ss << _cplex.getStatus();

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

