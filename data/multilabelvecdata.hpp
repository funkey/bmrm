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

#ifndef _MULTILABELVECDATA_HPP_
#define _MULTILABELVECDATA_HPP_

#include <vector>
#include <iostream>

#include "common.hpp"
#include "sml.hpp"
#include "data.hpp"
#include "vecfeature.hpp"
#include "varlenveclabel.hpp"


/** Container for dataset of sequence or variable length labels and vector features 
 * 
 *  The dataset consists of a list of examples (lines).  Each example
 *  consists of a scalar (or vector) valued label and a feature vector
 *  stored in different files. Also, special (e.g. query) id can be
 *  assigned to examples to form non-overlapping subsets of examples.
 *   
 *  See CVecFeature and CVarLenVecLabel for details.
 *
 * Notes:
 *   1.  The data loading starts during the object construction.
 */
class CMultilabelVecData : public CData, public CVecFeature, public CVarLenVecLabel
{        
   protected:
      /** To create matrix row view for feature matrices?
       */
      bool matrixRowView;
      
   public:
      /** Constructor.
       *  Divide data into "nparts" parts and load only the "start"-th part.
       */
      CMultilabelVecData(unsigned int start=0, unsigned int nparts=1);

      /** Destructor
       */
      virtual ~CMultilabelVecData() {}
      
      /** Whether data object has explicit access to each data point
       *  @return the value of matrixRowView
       */
      bool const HasMatrixRowView() 
      {
         return featureMatrixRowView;
      }
      
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
         return numOfExample; 
      }
      
      /** 
       * Return the total number of examples in this dataset.
       */
      virtual unsigned int size(void) const 
      { 
         return numOfAllExample; 
      }
      
      /** 
       * Return the dimension of example in this dataset.
       */
      virtual unsigned int dim(void) const 
      { 
         return featureDimension; 
      }
      
      /** 
       * Return true if we use a bias feature. False otherwise. 
       */
      virtual bool bias(void) const 
      { 
         return biasFlag; 
      }      
      
      virtual void ShowData()
      {
         ShowFeature();
         ShowLabel();
      }
};

#endif
