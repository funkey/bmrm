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

#ifndef _COSTFACTORY_HPP_
#define _COSTFACTORY_HPP_

#include "cost.hpp"

#include "hammingcost.hpp"
#include "zerocost.hpp"

/** Factory class to create cost functions. These functions can be used in some
 * loss functions.
 *
 * Cost functions give a measure of dissimilarity between two (usally vector-)
 * labes y and y'.
 */
class CostFactory {

public:

	static Cost* GetCost(CData* data) {

		Configuration &config = Configuration::GetInstance();

		if (config.IsSet("Loss.costFunctionType")) {

			std::string costFunctionType = config.GetString("Loss.costFunctionType");

			if (costFunctionType == "NORMALIZED_HAMMING") {

				// instantiate normalized Hamming cost function
				return new HammingCost(dynamic_cast<CVecLabel*>(data), true);

			} else if (costFunctionType == "HAMMING") {

				// instantiate Hamming cost function
				return new HammingCost(dynamic_cast<CVecLabel*>(data), false);

			} else if (costFunctionType == "ZERO") {

				// instantiate zero cost function
				return new ZeroCost();

			} else {

				throw CBMRMException(
						"ERROR: unrecognised cost function (" + costFunctionType + ")\n",
						"CostFactory::GetCost()");
			}
		} else
			throw CBMRMException(
					"ERROR: no cost function specified!\n",
					"CostFactory::GetCost()");
	}
};

#endif // _COSTFACTORY_HPP_

