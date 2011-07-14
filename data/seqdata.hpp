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

#ifndef _SEQDATA_HPP_
#define _SEQDATA_HPP_

#include <vector>
#include <iostream>

#include "common.hpp"
#include "sml.hpp"
#include "data.hpp"
#include "seqfeature.hpp"
#include "seqlabel.hpp"


/** Container for dataset of sequence labels and vectorial feature vectors
 * 
 *  Refer CSeqFeature and CSeqLabel for details.
 *
 * Notes:
 *   1.  The data loading starts during the object construction.
 *   2.  For centralized dataset and serial computation mode only!
 */
class CSeqData : public CData, public CSeqFeature, public CSeqLabel
{
public:
        CSeqData();
        virtual ~CSeqData() {}


        /** whether dataset contains labels
         */
        bool HasLabel() 
        { 
                return (Y.size() != 0); 
        }
        
        /** 
         * Return the number of examples in this sub-dataset.
         */
        virtual unsigned int slice_size(void) const 
        { 
                return numOfSeq; 
        }
      
        /** 
         * Return the total number of examples in this dataset.
         */
        virtual unsigned int size(void) const 
        { 
                return numOfAllSeq; 
        }
      
        /** 
         * Return the dimension of example in this dataset.
         */
        virtual unsigned int dim(void) const 
        { 
                return featureDimension; 
        }      
        
        /** No bias
         */
        virtual bool bias() const
        {
                return false;
        }
        
};

#endif
