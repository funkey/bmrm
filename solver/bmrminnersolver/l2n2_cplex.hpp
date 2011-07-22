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
 * Created: (07/20/2011)
 */

#ifndef _L2N2_CPLEX_HPP_
#define _L2N2_CPLEX_HPP_

#include <iostream>

#include <ilcplex/ilocplex.h>

#include "sml.hpp"
#include "model.hpp"
#include "l2n2_bmrmdualinnersolver.hpp"

class CL2N2_Cplex : public CL2N2_BMRMDualInnerSolver {

public:

	CL2N2_Cplex(double lambda);

protected:

	/** Solve QP
	*/
	virtual void SolveQP();

private:

	// the dimension x during the last call to SolveQP
	int _oldDim;

	// the Ilo environment
	IloEnv _env;

	// the solver variables x
	IloNumVarArray _vars;

	// the objective
	IloObjective _objective;

	// the expression of the objective
	IloNumExpr _expr;

	// the linear constraint
	IloRange _constraint;

	// the Ilo model containing the objective and constraints
	IloModel _model;

	// the Ilo Cplex solver object
	IloCplex _cplex;
};

#endif // _L2N2_CPLEX_HPP_

