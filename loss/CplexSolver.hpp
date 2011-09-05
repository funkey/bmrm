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

#ifndef _CPLEX_SOLVER_HPP_
#define _CPLEX_SOLVER_HPP_

#include <string>

#include "LinearProgramSolver.hpp"
#include <ilcplex/ilocplex.h>

using namespace std;

/**
 * Cplex interface to solve the following integer linear program:
 *
 * min  <a,x>
 * s.t. Ax  == b
 *      Cx  <= d
 *      x_i \in {0,1} for all i
 *
 * Where (A,b) describes all linear equality constraints, (C,d) all linear
 * inequality constraints and x is a binary vector. a is a real-valued vector
 * denoting the coefficients of the objective.
 */
class CplexSolver : public LinearProgramSolver {

public:

	CplexSolver(
		unsigned int numVariables,
		unsigned int numEqConstraints,
		unsigned int numIneqConstraints);

	virtual ~CplexSolver();

	/**
	 * Set the coefficients of the objective.
	 * @param a [read] A real valued vector of the coefficients.
	 */
	virtual void SetObjective(const TheMatrix& a, double constant, Sense sense);

	/**
	 * Set the linear equality constraints Ax == b.
	 * @param A [read] Coefficient matrix of the constraints.
	 * @param b [read] Values of the constraints.
	 */
	virtual void SetEqualities(const TheMatrix& A, const TheMatrix& b);

	/**
	 * Set the linear inequality constraints Ax <= b.
	 * @param A [read] Coefficient matrix of the constraints.
	 * @param b [read] Values of the constraints.
	 */
	virtual void SetInequalities(const TheMatrix& A, const TheMatrix& b);

	/**
	 * Solve the integer linear program.
	 * @param x     [write] The solution x of the integer linear program.
	 * @param value [write] The solution value of the objective.
	 * @param msg   [write] The return message of the solver.
	 * @return true if the program was solved successfully.
	 */
	virtual bool Solve(TheMatrix& x, double& value, string& msg);

private:

	/**
	 * Set the linear (in)equality constraints Ax <=|==|>= b.
	 * @param A [read] Coefficient matrix of the constraints.
	 * @param b [read] Values of the constraints.
	 * @param relation [read] The relation of the constraints:
	 *                          -1 --> "<="
	 *                           0 --> "=="
	 *                           1 --> ">="
	 */
	void SetConstraints(const TheMatrix& A, const TheMatrix& b, int relation);

	// the Ilo environment
	IloEnv _env;

	// the binary variables x
	IloNumVarArray _variables;

	// the objective
	IloObjective _objective;

	// the coefficients of the objective
	IloNumArray _coefs;

	// the Ilo model containing the objective and constraints
	IloModel _model;

	// the Ilo Cplex solver object
	IloCplex _cplex;

	// the verbosity of the output - set via Cplex.verbosity in config file
	int _verbosity;

	// a value by which to scale the objective
	double _scale;
};

#endif // _CPLEX_SOLVER_HPP_

