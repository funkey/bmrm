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
 * Created: (07/15/2011)
 */

#ifndef _VECCONSDATA_HPP_
#define _VECCONSDATA_HPP_

#include <vecdata.hpp>
#include <vecconstraints.hpp>

/** Container for dataset of vector features, vector labels and linear
 * constraints on the labels.
 * 
 *  The dataset consists of a list of examples (lines).  Each example
 *  consist of a scalar (or vector) valued label and a feature vector
 *  stored in different files. Additionally, linear constraints on the labels
 *  can be stored.
 */
class CConsVecData :
		public CVecData,
		public CVecConstraints {

public:

	CConsVecData(
			unsigned int start  = 0,
			unsigned int nparts = 1) {

		if (nparts != 1) {

			string msg =
					"CConsVecData does not support distribution of data (yet)";
			throw CBMRMException(msg, "CConsVecData::CConsVecData");
		}

		LoadConstraintData(NumOfLabel());
	}

	virtual ~CConsVecData() {}
};

#endif // _VECCONSDATA_HPP_

