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
	_vecconstraints_verbosity(0),
	_Aeq(0),
	_beq(0),
	_Aineq(0),
	_bineq(0),
	_numOfConstraints(0),
	_numOfEqualities(0),
	_numOfInequalities(0),
	_numOfVariables(0),
	_constraintsFile("")
{}

CVecConstraints::~CVecConstraints() {

	if (_Aeq)
		delete _Aeq;
	if (_beq)
		delete _beq;
	if (_Aineq)
		delete _Aineq;
	if (_bineq)
		delete _bineq;
}

void
CVecConstraints::LoadConstraintData(unsigned int numOfVariables) {

	_numOfVariables = numOfVariables;

	CTimer scanconstraintstime;
	CTimer loadconstraintstime;

	// get configurations
	Configuration &config = Configuration::GetInstance();

	_constraintsFile = config.GetString("Data.constraintsFile");
	if(config.IsSet("Data.verbosity"))
		_vecconstraints_verbosity = config.GetInt("Data.verbosity");

	// collect some properties of the dataset
	if(_vecconstraints_verbosity >= 1)
		cout << "Scanning constraints file... "<< endl;

	scanconstraintstime.Start();
	ScanConstraintsFile();
	scanconstraintstime.Stop();

	// read constraintss into memory
	if(_vecconstraints_verbosity >= 1)
		cout << "Loading constraints file... "<< endl;

	loadconstraintstime.Start();
	LoadConstraints();
	loadconstraintstime.Stop();

	if(_vecconstraints_verbosity >= 2) {
		cout << "scanconstraintstime	: " << scanconstraintstime.CPUTotal() << endl;
		cout << "loadconstraintstime	: " << loadconstraintstime.CPUTotal() << endl;
	}
}

void
CVecConstraints::ScanConstraintsFile() {

	string line = "";
	string token = "";
	ifstream constraintsFp;

	_numOfConstraints  = 0;
	_numOfEqualities   = 0;
	_numOfInequalities = 0;

	unsigned int numOfVariables    = 0;

	constraintsFp.open(_constraintsFile.c_str());
	if(!constraintsFp.good()) {

		string msg = "Cannot open constraints file <" + _constraintsFile + ">!";
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

				_numOfInequalities++;
				break;
			}

			if (token[0] == '=') {

				_numOfEqualities++;
				break;
			}

			double coefficient;
			unsigned int id;

			ScanExpression(token, coefficient, id);

			numOfVariables = max(numOfVariables, id+1);
		}
	}

	if (numOfVariables > _numOfVariables) {

		string msg =
				string("Number of variables in constraint file (") +
				lexical_cast<string>(numOfVariables) +
				") is more than expected (" +
				lexical_cast<string>(_numOfVariables);
		throw CBMRMException(msg, "CVecConstraints::ScanConstraintsFile()");
	}

	_numOfConstraints = _numOfEqualities + _numOfInequalities;

	if (_vecconstraints_verbosity >= 1)
		if(_numOfConstraints == 0)
			cout << "Constraints file is empty" << endl;
		else
			cout << _numOfConstraints << " constraints found ("
			     << _numOfEqualities << " eq., "
			     << _numOfInequalities << " ineq.)" << endl;

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
		             _constraintsFile + ">. Expected *, got \"" +
		             op + "\".";
		throw CBMRMException(msg, "CVecConstraints::ScanExpression()");
	}
}

void
CVecConstraints::LoadConstraints() {

	if (_numOfEqualities > 0) {

		_Aeq   = new TheMatrix(_numOfEqualities, _numOfVariables, SML::SPARSE);
		_beq   = new TheMatrix(_numOfEqualities, 1, SML::DENSE);
	}

	if (_numOfInequalities > 0) {

		_Aineq = new TheMatrix(_numOfInequalities, _numOfVariables, SML::SPARSE);
		_bineq = new TheMatrix(_numOfInequalities, 1, SML::DENSE);
	}

	string line = "";
	string token = "";
	ifstream constraintsFp;

	constraintsFp.open(_constraintsFile.c_str());
	if(!constraintsFp.good()) {

		string msg = "Cannot open constraints file <" + _constraintsFile + ">!";
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
					_Aineq->Set(nextInequality, ids[i], coefs[i]);

				double value;
				iss >> value;

				_bineq->Set(nextInequality, value);

				nextInequality++;
				break;
			}

			if (token[0] == '>') {

				for (unsigned int i = 0; i < ids.size(); i++)
					_Aineq->Set(nextInequality, ids[i], -coefs[i]);

				double value;
				iss >> value;

				_bineq->Set(nextInequality, -value);

				nextInequality++;
				break;
			}

			if (token[0] == '=') {

				for (unsigned int i = 0; i < ids.size(); i++)
					_Aeq->Set(nextEquality, ids[i], coefs[i]);

				double value;
				iss >> value;

				_beq->Set(nextEquality, value);

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

	if (nextEquality != _numOfEqualities) {

		string msg =
				string("Only ") + lexical_cast<string>(nextEquality) +
				" equalities read, " + lexical_cast<string>(_numOfEqualities) +
				" expected.";

		throw CBMRMException(msg, "CVecConstraints::LoadConstraints()");
	}
	if (nextInequality != _numOfInequalities) {

		string msg =
				string("Only ") + lexical_cast<string>(nextInequality) +
				" inequalities read, " + lexical_cast<string>(_numOfInequalities) +
				" expected.";

		throw CBMRMException(msg, "CVecConstraints::LoadConstraints()");
	}

	if (_vecconstraints_verbosity >= 2) {

		if (_numOfEqualities > 0) {

			cout << "Equalities: " << endl;
			_Aeq->Print();
			_beq->Print();
		}

		if (_numOfInequalities > 0) {

			cout << "Inequalities: " << endl;
			_Aineq->Print();
			_bineq->Print();
		}
	}

	constraintsFp.close();
}
