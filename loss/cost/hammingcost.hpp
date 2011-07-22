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
 * Created: (07/21/2011)
 */

#ifndef _HAMMINGCOST_HPP_
#define _HAMMINGCOST_HPP_

#include "cost.hpp"
#include "veclabel.hpp"

/** The Hamming cost function. The default implementation gives
 *
 *   Δ(y,y') = 1/N \sum_i | y_i - y'_i |
 *
 * If the "normalize" flag is given, N will be the the maximal Hamming cost,
 * otherwise it is just 1.
 *
 * Via the configuration options HammingCost.zeroToOneCost and
 * HammingCost.oneToZeroCost the Hamming cost can be made asymmetric.
 */
class HammingCost : public Cost {

public:

	/** Default constructor.
	 *
	 * @param data      [read] Pointer to the label data.
	 * @param normalize [read] Whether to normalize the Hamming distance.
	 */
	HammingCost(
			CVecLabel* data,
			bool normalize);

	/** Get the constant contribution of the cost function to the objective.
	 *
	 * @param y [read]  The ground truth label.
	 * @param c [write] The constant term of the cost function.
	 */
	virtual void constantContribution(const TheMatrix& y, double& c) const;

	/** Get the linear contribution of the cost function to the objective.
	 *
	 * @param y [read]  The ground truth label.
	 * @param a [write] The linear coefficients.
	 */
	virtual void linearContribution(const TheMatrix& y, TheMatrix& a) const;

private:

	// the number of variables in a label y
	int _numVariables;

	// normalize the Hamming cost?
	bool _normalize;

	// deviation costs from reference label y
	double _oneToZeroCost;
	double _zeroToOneCost;
};

#endif // _HAMMINGCOST_HPP_

