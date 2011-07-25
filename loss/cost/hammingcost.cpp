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

#include "hammingcost.hpp"
#include "configuration.hpp"

HammingCost::HammingCost(CVecLabel* data, bool normalize) :
	_normalize(normalize),
	_oneToZeroCost(1.0),
	_zeroToOneCost(1.0) {

	// if the labels have dimension 1, we assume that the Hamming distance is
	// supposed to be computed over all labels (everything else would be
	// pointless)
	if (data->LabelDimension() == 1)
		_numVariables = data->NumOfLabel();
	else
		_numVariables = data->LabelDimension();

	// get the configuration
	Configuration &config = Configuration::GetInstance();

	if (config.IsSet("HammingCost.oneToZeroCost"))
		_oneToZeroCost = config.GetDouble("HammingCost.oneToZeroCost");
	if (config.IsSet("HammingCost.zeroToOneCost"))
		_zeroToOneCost = config.GetDouble("HammingCost.zeroToOneCost");

}

void
HammingCost::constantContribution(const TheMatrix& y, double& c) const {

	int numOnes = 0;
	c = 0;

	for (int i = 0; i < _numVariables; i++) {

		double value;

		// accumulate constant cost term
		y.Get(i, value);
		numOnes += value;

		c += _oneToZeroCost*value;
	}

	// normalize Hamming distance?
	if (_normalize) {

		int numZeros = _numVariables - numOnes;

		// maximum Hamming distance
		double maxCost =
				numOnes*_oneToZeroCost +
				numZeros*_zeroToOneCost;

		c *= 1.0/maxCost;
	}
}

void
HammingCost::linearContribution(const TheMatrix& y, TheMatrix& a) const {

	int numOnes = 0;

	for (int i = 0; i < _numVariables; i++) {

		double value;

		// calculate linear cost coefficients
		y.Get(i, value);
		a.Set(i, (_zeroToOneCost - (_zeroToOneCost + _oneToZeroCost)*value));

		numOnes += value;
	}

	// normalize Hamming distance?
	if (_normalize) {

		int numZeros = _numVariables - numOnes;

		// maximum Hamming distance
		double maxCost =
				numOnes*_oneToZeroCost +
				numZeros*_zeroToOneCost;

		a.Scale(1.0/maxCost);
	}
}
