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
#include <fstream>
#include <vecconstraints.hpp>
#include <boost/lexical_cast.hpp>

#include "common.hpp"
#include "configuration.hpp"
#include "timer.hpp"

using namespace std;
using namespace boost;

CVecConstraints::CVecConstraints() :
	vecconstraints_verbosity(0),
	Aeq(0),
	beq(0),
	Aineq(0),
	bineq(0),
	numOfConstraints(0),
	numOfEqualities(0),
	numOfInequalities(0),
	numOfVariables(0),
	constraintsFile("")
{}

CVecConstraints::~CVecConstraints() {

	if (Aeq)
		delete Aeq;
	if (beq)
		delete beq;
	if (Aineq)
		delete Aineq;
	if (bineq)
		delete bineq;
}

void
CVecConstraints::LoadConstraintData() {

	CTimer scanconstraintstime;
	CTimer loadconstraintstime;

	// get configurations
	Configuration &config = Configuration::GetInstance();

	constraintsFile = config.GetString("Data.constraintsFile");
	if(config.IsSet("Data.verbosity"))
		vecconstraints_verbosity = config.GetInt("Data.verbosity");

	// collect some properties of the dataset
	if(vecconstraints_verbosity >= 1)
		cout << "Scanning constraints file... "<< endl;

	scanconstraintstime.Start();
	ScanConstraintsFile();
	scanconstraintstime.Stop();

	// read constraintss into memory
	if(vecconstraints_verbosity >= 1)
		cout << "Loading constraints file... "<< endl;

	loadconstraintstime.Start();
	LoadConstraints();
	loadconstraintstime.Stop();

	if(vecconstraints_verbosity >= 2) {
		cout << "scanconstraintstime	: " << scanconstraintstime.CPUTotal() << endl;
		cout << "loadconstraintstime	: " << loadconstraintstime.CPUTotal() << endl;
	}
}

void
CVecConstraints::ScanConstraintsFile() {

	string line = "";
	string token = "";
	ifstream constraintsFp;

	numOfConstraints  = 0;
	numOfEqualities   = 0;
	numOfInequalities = 0;
	numOfVariables    = 0;

	constraintsFp.open(constraintsFile.c_str());
	if(!constraintsFp.good()) {

		string msg = "Cannot open constraints file <" + constraintsFile + ">!";
		throw CBMRMException(msg, "CVecConstraints::ScanConstraintsFile()");
	}

	while(!constraintsFp.eof()) {

		getline(constraintsFp, line);

		trim(line);
		if(IsBlankLine(line)) continue;  // blank line
		if(line[0] == '#') continue;  // comment line

		istringstream iss(line);

		while (!iss.eof()) {

			iss >> token;

			if (token[0] == '<' || token[0] == '>') {

				numOfInequalities++;
				break;
			}

			if (token[0] == '=') {

				numOfEqualities++;
				break;
			}

			double coefficient;
			unsigned int id;

			ScanExpression(token, coefficient, id);

			numOfVariables = max(numOfVariables, id+1);
		}
	}

	numOfConstraints = numOfEqualities + numOfInequalities;

	if (vecconstraints_verbosity >= 1)
		if(numOfConstraints == 0)
			cout << "Constraints file is empty" << endl;
		else
			cout << numOfConstraints << " constraints found ("
			     << numOfEqualities << " eq., "
			     << numOfInequalities << " ineq.)" << endl;

	constraintsFp.close();
}

void
CVecConstraints::ScanExpression(string expr, double& coefficient, unsigned int& id) {

	istringstream iss(expr);

	char op;

	iss >> coefficient;
	iss >> op;
	iss >> id;

	if (op != '*') {
		string msg = "Malformed constraints file <" +
		             constraintsFile + ">. Expected *, got \"" +
		             op + "\".";
		throw CBMRMException(msg, "CVecConstraints::ScanExpression()");
	}
}

void
CVecConstraints::LoadConstraints() {

	Aeq   = new TheMatrix(numOfEqualities, numOfVariables, SML::SPARSE);
	beq   = new TheMatrix(numOfEqualities, 1, SML::DENSE);
	Aineq = new TheMatrix(numOfInequalities, numOfVariables, SML::SPARSE);
	bineq = new TheMatrix(numOfInequalities, 1, SML::DENSE);

	string line = "";
	string token = "";
	ifstream constraintsFp;

	constraintsFp.open(constraintsFile.c_str());
	if(!constraintsFp.good()) {

		string msg = "Cannot open constraints file <" + constraintsFile + ">!";
		throw CBMRMException(msg, "CVecConstraints::LoadConstraints()");
	}

	unsigned int nextEquality   = 0;
	unsigned int nextInequality = 0;

	while(!constraintsFp.eof()) {

		getline(constraintsFp, line);

		trim(line);
		if(IsBlankLine(line)) continue;  // blank line
		if(line[0] == '#') continue;  // comment line

		istringstream iss(line);

		vector<double>       coefs;
		vector<unsigned int> ids;

		while (!iss.eof()) {

			iss >> token;

			if (token[0] == '<') {

				for (unsigned int i = 0; i < ids.size(); i++)
					Aineq->Set(nextInequality, ids[i], coefs[i]);

				double value;
				iss >> value;

				bineq->Set(nextInequality, value);

				nextInequality++;
				break;
			}

			if (token[0] == '>') {

				for (unsigned int i = 0; i < ids.size(); i++)
					Aineq->Set(nextInequality, ids[i], -coefs[i]);

				double value;
				iss >> value;

				bineq->Set(nextInequality, -value);

				nextInequality++;
				break;
			}

			if (token[0] == '=') {

				for (unsigned int i = 0; i < ids.size(); i++)
					Aeq->Set(nextEquality, ids[i], coefs[i]);

				double value;
				iss >> value;

				beq->Set(nextEquality, value);

				nextEquality++;
				break;
			}

			double coefficient;
			unsigned int id;

			ScanExpression(token, coefficient, id);

			coefs.push_back(coefficient);
			ids.push_back(id);
		}
	}

	if (nextEquality != numOfEqualities) {

		string msg =
				string("Only ") + lexical_cast<string>(nextEquality) +
				" equalities read, " + lexical_cast<string>(numOfEqualities) +
				" expected.";

		throw CBMRMException(msg, "CVecConstraints::LoadConstraints()");
	}
	if (nextInequality != numOfInequalities) {

		string msg =
				string("Only ") + lexical_cast<string>(nextInequality) +
				" inequalities read, " + lexical_cast<string>(numOfInequalities) +
				" expected.";

		throw CBMRMException(msg, "CVecConstraints::LoadConstraints()");
	}

	if (vecconstraints_verbosity >= 2) {

		cout << "Equalities: " << endl;
		Aeq->Print();
		beq->Print();
		cout << "Inequalities: " << endl;
		Aineq->Print();
		bineq->Print();
	}

	constraintsFp.close();
}
