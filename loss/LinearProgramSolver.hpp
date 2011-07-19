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

#ifndef _LINEAR_PROGRAM_SOLVER_HPP_
#define _LINEAR_PROGRAM_SOLVER_HPP_

#include <string>

#include <sml.hpp>

using namespace std;

/**
 * Abstract class for linear program solvers. Implementations are supposed to
 * solve the following integer linear program:
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
class LinearProgramSolver {

public:

	/** Used to indicate whether the solver should attempt to minimize or
	 * maximize the objective.
	 */
	enum Sense {

		MINIMIZE,
		MAXIMIZE
	};

	LinearProgramSolver(
		unsigned int numVariables,
		unsigned int numEqConstraints,
		unsigned int numIneqConstraints) :
		_numVariables(numVariables),
		_numEqConstraints(numEqConstraints),
		_numIneqConstraints(numIneqConstraints)
		{}

	virtual ~LinearProgramSolver() {}

	/**
	 * Set the coefficients of the objective.
	 * @param a [read] A real valued vector of the coefficients.
	 */
	virtual void SetObjective(
			const TheMatrix& a,
			double constant,
			Sense sense) = 0;

	/**
	 * Set the linear equality constraints Ax == b.
	 * @param A [read] Coefficient matrix of the constraints.
	 * @param b [read] Values of the constraints.
	 */
	virtual void SetEqualities(const TheMatrix& A, const TheMatrix& b) = 0;

	/**
	 * Set the linear inequality constraints Cx <= d.
	 * @param C [read] Coefficient matrix of the constraints.
	 * @param d [read] Values of the constraints.
	 */
	virtual void SetInequalities(const TheMatrix& C, const TheMatrix& d) = 0;

	/**
	 * Solve the integer linear program.
	 * @param x     [write] The solution x of the integer linear program.
	 * @param value [write] The solution value of the objective.
	 * @param msg   [write] The return message of the solver.
	 * @return true if the program was solved successfully.
	 */
	virtual bool Solve(TheMatrix& x, double& value, string& msg) = 0;

protected:

	unsigned int _numVariables;

	unsigned int _numEqConstraints;

	unsigned int _numIneqConstraints;

};

#endif // _LINEAR_PROGRAM_SOLVER_HPP_

