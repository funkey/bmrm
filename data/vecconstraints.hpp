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

#ifndef _VECCONSTRAINTS_HPP_
#define _VECCONSTRAINTS_HPP_

#include <vector>
#include <iostream>
#include <string>

#include "sml.hpp"

using namespace std;

/** Container for linear constraints on the structured labels.
 *
 *   Label file format: rows of <constraint>
 *   where
 *   <constraint> .=. <expr> <rel> <value>
 *   <expr>       .=. <coef_1>*<id_1> [ <coef_2>*<id_2> ... <coef_k>*<id_k> ]
 *   <rel>        .=. <= | == | =>
 *   <value>      .=. scalar right hand side value of constraint
 *   <coef_i>     .=. coefficient for variable id_i
 *   <id_i>       .=. number of the variable as given in label and feature
 *                    files
 *
 * Notes:
 *   1.  The data loading starts during the object construction.
 *   2.  For centralized dataset and serial computation mode only!
 */
class CVecConstraints {

public:

	/** Constructor
	 */
	CVecConstraints();

	/** Destructor
	 */
	virtual ~CVecConstraints();

	/**
	 * Scans and loads the constraints file.
	 * @param numVariables [read] The expected number of variables involved in
	 *                            the constraints. This will influence the shape
	 *                            of the resulting coefficient matrices.
	 */
	virtual void LoadConstraintData(unsigned int numVariables);

	bool HasData() {

		return (_Aeq != 0);
	}

	unsigned int NumOfConstraints() {

		return _numOfConstraints;
	}

	unsigned int NumOfEqualities() {

		return _numOfEqualities;
	}

	unsigned int NumOfInequalities() {

		return _numOfInequalities;
	}

	const TheMatrix& GetEqualityCoefs() {

		return *_Aeq;
	}

	const TheMatrix& GetEqualityValues() {

		return *_beq;
	}

	const TheMatrix& GetInequalityCoefs() {

		return *_Aineq;
	}

	const TheMatrix& GetInequalityValues() {

		return *_bineq;
	}

	std::string ConstraintsFile() {

		return _constraintsFile;
	}

	void ShowConstraints() {

		_Aeq->Print();
		_beq->Print();
		_Aineq->Print();
		_bineq->Print();
	}

protected:

	/** Scan an expression of the form <coef>*<id>
	 */
	void ScanExpression(string expr, double& coefficient, unsigned int& id);

	/** Scan constraints file to determine number of (in)equalities
	 */
	void ScanConstraintsFile();

	/** Allocate constraints matrices and vectors and load constraints from file
	 */
	void LoadConstraints();

	/** Verbosity level
	 */
	int _vecconstraints_verbosity;

	/** Coefficient matrix of equalities (each row is a constraint)
	 */
	TheMatrix* _Aeq;

	/** Values vector of equalities
	 */
	TheMatrix* _beq;

	/** Coefficient matrix of inequalities (each row is a constraint)
	 */
	TheMatrix* _Aineq;

	/** Values vector of inequalities
	 */
	TheMatrix* _bineq;

	/** Total number of constraints.
	 */
	unsigned int _numOfConstraints;

	/** Number of equalities.
	 */
	unsigned int _numOfEqualities;

	/** Number of inequalities.
	 */
	unsigned int _numOfInequalities;

	/** The number of variables involved in constraints.
	 */
	unsigned int _numOfVariables;

	/** Name of the file containing constraints
	 */
	std::string _constraintsFile;
};

#endif // _VECCONSTRAINTS_HPP_

