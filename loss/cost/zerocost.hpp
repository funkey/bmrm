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

#ifndef _ZEROCOST_HPP_
#define _ZEROCOST_HPP_

#include "cost.hpp"
#include "veclabel.hpp"

/** A zero cost function.
 *
 *   Δ(y,y') = 0
 */
class ZeroCost : public Cost {

public:

	/** Default constructor.
	 */
	ZeroCost() {};

	/** Get the constant contribution of the cost function to the objective.
	 *
	 * @param y [read]  The ground truth label.
	 * @param c [write] The constant term of the cost function.
	 */
	virtual void constantContribution(const TheMatrix& y, double& c) const {
	
		c = 0.0;
	};

	/** Get the linear contribution of the cost function to the objective.
	 *
	 * @param y [read]  The ground truth label.
	 * @param a [write] The linear coefficients.
	 */
	virtual void linearContribution(const TheMatrix& y, TheMatrix& a) const {
	
		a.Zero();
	}
};

#endif // _ZEROCOST_HPP_

