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
 * Created: (30/01/2008) 
 *
 * Last Updated:
 */

#ifndef _SEQDATA_CPP_
#define _SEQDATA_CPP_

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
#include "seqdata.hpp"

   
/**  Constructor
 */
CSeqData::CSeqData()
   : CData(), 
     CSeqFeature(), 
     CSeqLabel()
{
        // sanity check
        if(numOfAllSeq != numOfAllLabel)
        {
                std::ostringstream msg;
                msg << "Number of sequences (" << numOfAllSeq << ") does not match number of labels (" << numOfAllLabel << ")";
                throw CBMRMException(msg.str(),"CSeqData::CSeqData()");
        }
        
        // dataset statistics
        if(verbosity)
        {
                std::cout << "Dataset properties:"  << std::endl;
                std::cout << "1.  Feature file              : " << featureFile << std::endl;
                std::cout << "2.  Label file                : " << labelFile << std::endl;
                std::cout << "3.  Number of sequences       : " << numOfAllSeq << std::endl;    
                std::cout << "4.  Number of nonzero features: " << nnz << std::endl;
                std::cout << "5.  Max sequence length       : " << maxSeqLen << std::endl;
                std::cout << "6.  Min sequence length       : " << minSeqLen << std::endl;
                std::cout << "7.  Total sequence length     : " << totalSeqLen << std::endl;
                std::cout << "8.  Avg. sequence length      : " << (double)totalSeqLen/numOfAllSeq << std::endl;
                std::cout << "9.  Dim. of feature vector    : " << featureDimension << std::endl; 
                //std::cout << "8.  Dataset density           : " << density*100.0 << std::endl;
                //std::cout << "9.  Avg. nnz(feature vec.)    : " << ((double)nnz)/numOfAllSeq << std::endl;
                std::cout << std::endl;
        }
}

#endif
