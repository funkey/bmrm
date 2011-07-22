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

#include <iostream>
#include <sstream>

#include <ilcplex/ilocplex.h>

#include "l2n2_cplex.hpp"

using namespace std;

CL2N2_Cplex::CL2N2_Cplex(double lambda) :
	CL2N2_BMRMDualInnerSolver(lambda) {}

void
CL2N2_Cplex::SolveQP() {

	// link to the paper:
	//
	// Q = 1/λA^T*A
	// f = -b
	// l … lower bounds on alpha
	// u … lower bounds on alpha
	// a … constraint matrix (one row of 1s)
	// b … right hand side of constraints (=1)
	//
	// -> we have to minimize g = 0.5*x^T*Q*x + x^T*f

	if (verbosity > 2) {

		cout << "[L2N2_Cplex::SolveQP] Q:" << endl;

		for (int i = 0; i < dim; i++) {

			for (int j = 0; j < dim; j++) {

				cout << Q[j + i*dim] << " ";
			}

			cout << endl;
		}

		cout << "[L2N2_Cplex::SolveQP] f:" << endl;

		for (int i = 0; i < dim; i++)
			cout << f[i] << " ";
		cout << endl;

		cout << "[L2N2_Cplex::SolveQP] a:" << endl;

		for (int i = 0; i < dim; i++)
			cout << a[i] << " ";
		cout << endl;

		cout << "[L2N2_Cplex::SolveQP] b:" << endl << *b << endl;

		cout << "[L2N2_Cplex::SolveQP] l:" << endl;

		for (int i = 0; i < dim; i++)
			cout << l[i] << " ";
		cout << endl;

		cout << "[L2N2_Cplex::SolveQP] u:" << endl;

		for (int i = 0; i < dim; i++)
			cout << u[i] << " ";
		cout << endl;
	}

	if (verbosity > 1)
		cout << "[L2N2_Cplex::SolveQP] allocating variable memory" << endl;

	IloEnv env;

	IloModel model(env);

	IloNumVarArray vars(env, dim, 0.0, 1.0, ILOFLOAT);

	if (verbosity > 1)
		cout << "[L2N2_Cplex::SolveQP] " << dim << " variables allocated" << endl;

	/////////////////////
	// setup objective //
	/////////////////////

	if (verbosity > 1)
		cout << "[L2N2_Cplex::SolveQP] setting objective" << endl;

	IloObjective objective(env);
	objective.setSense(IloObjective::Minimize);

	IloNumExpr expr(env);

	try {

		// quadratic term 0.5*x^T*Q*x
		for (int row = 0; row < dim; row++) {

			for (int col = row; col < dim; col ++) {

				double value = Q[col + row*dim];

				if (value == 0)
					continue;

				if (col == row) 
					expr = expr + static_cast<IloNum>(0.5*value)*vars[row]*vars[row];
				else
					expr = expr + static_cast<IloNum>(value)*vars[row]*vars[col];
			}
		}

		// linear term
		for (int i = 0; i < dim; i++)
			expr = expr + static_cast<IloNum>(f[i])*vars[i];

	} catch (IloMemoryException e) {

		cerr << "[L2N2_Cplex::SolveQP] exception: " << e.getMessage() << endl;
		throw CBMRMException(e.getMessage(), "L2N2_Cplex::SolveQP");
	}

	objective.setExpr(expr);

	model.add(objective);

	if (verbosity > 1)
		cout << "[L2N2_Cplex::SolveQP] objective set" << endl;

	//////////////////////
	// setup constraint //
	//////////////////////

	if (verbosity > 1)
		cout << "[L2N2_Cplex::SolveQP] setting constraint" << endl;

	IloRange constraint(env, *b, *b);

	for (int i = 0; i < dim; i++)
		constraint.setLinearCoef(vars[i], static_cast<IloNum>(a[i]));

	model.add(constraint);

	if (verbosity > 1)
		cout << "[L2N2_Cplex::SolveQP] constraint set" << endl;

	if (verbosity > 2)
		cout << "[L2N2_Cplex::SolveQP] Cplex model is:" << endl << model << endl;

	//////////////////
	// solve the QP //
	//////////////////

	if (verbosity > 1)
		cout << "[L2N2_Cplex::SolveQP] solving problem" << endl;

	IloCplex cplex(model);

	cplex.solve();

	// get solver result message
	stringstream ss;
	ss << cplex.getStatus() << flush;

	if (cplex.getStatus() != IloAlgorithm::Optimal)
		throw new CBMRMException(ss.str(), "L2N2_Cplex::SolveQP");

	if (verbosity > 0)
		cout << "[L2N2_Cplex::SolveQP] solver finished successfully" << endl;

	// extract solution
	IloNumArray values(env);

	cplex.getValues(values, vars);

	for (int i = 0; i < dim; i++)
		x[i] = static_cast<double>(values[i]);

	if (verbosity > 2) {

		cout << "[L2N2_Cplex::SolveQP] x: " << endl;

		for (int i = 0; i < dim; i++)
			cout << x[i] << " ";
		cout << endl;
	}

	// clean up
	values.end();
	cplex.end();
	constraint.end();
	expr.end();
	objective.end();
	vars.end();
	model.end();
	env.end();
}
