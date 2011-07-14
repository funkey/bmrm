/* Copyright (c) 2009, NICTA
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
 * Authors: Choon Hui Teo (ChoonHui.Teo@anu.edu.au)
 *
 * Created: (26/01/2008) 
 *
 * Last Updated:
 */

#ifndef _MULTILABELVECDATA_CPP_
#define _MULTILABELVECDATA_CPP_

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "common.hpp"
#include "sml.hpp"
#include "configuration.hpp"
#include "bmrmexception.hpp"
#include "timer.hpp"
#include "multilabelvecdata.hpp"

using namespace std;
   
/**  Constructor
 */
CMultilabelVecData::CMultilabelVecData(unsigned int start, unsigned int nparts)
   : CData(), 
     CVecFeature(), 
     CVarLenVecLabel()
{
   LoadFeatureMatrix(start, nparts);
   LoadVarLenLabel(startExampleIndex, numOfExample);
   
   // sanity check
   if(numOfAllExample != numOfAllLabel)
   {
      ostringstream msg;
      msg << "Number of examples (" << numOfAllExample << ") does not match number of labels (" << numOfAllLabel << ")";
      throw CBMRMException(msg.str(),"CMultilabelVecData::CMultilabelVecData()");
   }
   
   // dataset statistics
   if(verbosity)
   {
      cout << "Dataset properties:"  << endl;
      cout << "1.  Feature file              : " << featureFile << endl;
      cout << "2.  Label file                : " << labelFile << endl;
      cout << "3.  Number of subsets         : " << numOfAllSubset << endl;
      cout << "4.  Number of examples        : " << numOfAllExample << endl;    
      cout << "5.  Number of nonzero features: " << nnz; 
      if(biasFlag)
         cout << "(+ no. of examples for shifted hyperplane)";
      cout << endl;
      cout << "7.  Average label dimension   : " << avgLabelDimension << endl;
      cout << "8.  Dim. of feature vector    : " << featureDimension;
      if(biasFlag)
         cout << "(+ 1 for shifted hyperplane)";
      cout << endl; 
      cout << "9.  Dataset density           : " <<  density*100.0 << endl;
      cout << "10. Avg. nnz(feature vec.)    : " << ((double)nnz)/numOfAllExample << endl;
      cout << endl;
   }
}

#endif
