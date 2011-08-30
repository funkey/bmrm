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

#ifndef _SOFTMARGINLOSS_HPP_
#define _SOFTMARGINLOSS_HPP_

#include "sml.hpp"
#include "data.hpp"
#include "vecconsdata.hpp"
#include "loss.hpp"
#include "model.hpp"
#include "CplexSolver.hpp"
#include "cost.hpp"

class SoftMarginLoss : public CLoss {

public:

	/** Default constructor.
	 */
	SoftMarginLoss(CModel* model, CConsVecData* data);

	/** Destructor.
	 */
	~SoftMarginLoss();

	/** Compute the loss value given the weight vector w in model object  
	 * 
	 *  @param loss [write] value of the loss.
	 */
	void ComputeLoss(double& loss);


	/** Compute the loss value as well as the gradient of the loss function w.r.t weight 
	 *  vector in model object.
	 * 
	 *  @param loss [write] value of the loss.
	 *  @param grad [write] value of the gradient.
	 */
	void ComputeLossAndGradient(double& loss, TheMatrix& grad);


	/** Perform performance evaluation (with test labels) given the model object.
	 * 
	 *  @param model [read] model object containing weight vector (and possibly other extra parameter)
	 */
	void Evaluate(CModel *model) {
		// empty
	}

private:

	/** Computes the constant (c) and linear (l) contribution of the
	 * missclassification cost with respect to a ground-truth labeling.
	 * The resulting costs for any labeling y is then c + <l,y>.
	 */
	void ComputeCostContribution(const TheMatrix& groundTruth);

	/** Computes the constant (c) and linear (l) contribution of the
	 * gamma function with respect to a ground-truth labeling.
	 * The resulting gamma function for any labeling y is then c + <l,y>.
	 */
	void ComputeGammaContribution(const TheMatrix& groundTruth);

	/** Computes the number of auxilary variables and corresponding constraints
	 * needed to express a linear gamma function in the objective. Linear gamma
	 * functions result in quadratic terms and therefore auxilary variables are
	 * needed to restore a linear program.
	 */
	void AddAuxiliaryVariables();

	/** Determines and sets all linear constraints for the objective. Linear
	 * constraints may come from the data (explicitly given by the user) and/or
	 * from the auxilary variables that are used to express quadratic terms in
	 * the objective.
	 */
	void SetLinearConstraints();

	/** Checks whether the ground-truth fulfills the requirements imposed by the
	 * linear constraints.
	 */
	void CheckSolutionIntegrity(
			const TheMatrix& A,
			const TheMatrix& y,
			int relation,
			const TheMatrix& b);

	/** The training data.
	 */
	CConsVecData* _data;

	/** The model to train, i.e., the weight vector.
	 */
	CModel* _model;

	/** The number of components of the feature vectors in _data.
	 */
	unsigned int _numFeatures;

	/** The number of all binary variables in the objective.
	 */
	unsigned int _numVariables;

	/** The number of auxilary variables (used for modelling quadratic terms in
	 * the objective) in the objective.
	 */
	unsigned int _numAuxiliaryVariables;

	/** The number of all linear constraints on the objective.
	 */
	unsigned int _numEqualities;
	unsigned int _numInequalities;

	/** The number of linear constraints on the auxilary variables only.
	 */
	unsigned int _numAuxiliaryEqualities;
	unsigned int _numAuxiliaryInequalities;

	/** The linear solver used to compute the gradient and loss.
	 */
	CplexSolver* _solver;

	/** The cost function Δ(y,y') to be used for the margin scaling.
	 */
	Cost* _costFunction;

	/** The coefficient of the cost term.
	 */
	double _costFactor;

	/** The components of the loss function.
	 */
	double    _g_c;
	TheMatrix _g_l;
	double    _m_c;
	TheMatrix _m_l;
	double    _c_c;
	TheMatrix _c_l;

	/** Is the gamma function constant?
	 */
	bool _gammaConst;

	/** The ground-truth labeling as given by the data.
	 */
	TheMatrix _y;

	/** The verbosity of the output. Set via Loss.verbosity in the config file.
	 */
	int _verbosity;
};

#endif // _SOFTMARGINLOSS_HPP_

