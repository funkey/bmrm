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
 * Authors: Qinfeng Shi (qinfeng.shi@anu.edu.au)
 *          Choon Hui Teo (ChoonHui.Teo@anu.edu.au)
 *
 * Created: (14/04/2008) 
 *
 * Last Updated: (28/10/2008)
 */

#ifndef _SEQMULTICLASSDATA_HPP_
#define _SEQMULTICLASSDATA_HPP_

#include <vector>
#include <iostream>

#include "common.hpp"
#include "sml.hpp"
#include "data.hpp"
#include "seqmulticlassfeature.hpp"
#include "seqmulticlasslabel.hpp"


/** Container for dataset of sequence labels and vectorial feature vectors
 * 
 *  Refer CSeqFeature and CSeqLabel for details.
 *
 * Notes:
 *   1.  The data loading starts during the object construction.
 *   2.  For centralized dataset and serial computation mode only!
 */
namespace SMM {
    enum MARK{TRAIN_DATA,TEST_DATA,VALID_DATA};
}

class CSeqMulticlassData : public CData, public CSeqMulticlassFeature, public CSeqMulticlassLabel
{
	
   public:

      std::vector<int> cvmark;

      CSeqMulticlassData();
      virtual ~CSeqMulticlassData() {}
      
      /** set n-fold cross validation (test, train, validation sets marks) on examples 		
       */
      virtual void SetCrossValidationData(unsigned int numOfFold, unsigned int foldIndex);
      
      /** Extend the Raw Feature dimension to Tensor Feature dimension 
       */
      virtual void ExtendFeatures();
      
      /** Tensor the raw Phi1 feature to Phi1
       */
      virtual void TensorPhi1(TheMatrix* phi1, unsigned int classID, unsigned int personID, TheMatrix *v);
      
      /** Tensor the Phi2 feature to Phi2, Phi3
       */
      virtual void TensorPhi2(const TheMatrix*  phi2, unsigned int classIDPrev, unsigned int classID, unsigned int personID, int verbosity, TheMatrix *v);
      
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
       * Return the dimension of raw feature example in this dataset.
       */
      virtual unsigned int dim(void) const 
      { 
         return featureDimension; 
      }
      
      /** 
       * Return the dimension of tensor feature example in this dataset.
       */
      virtual unsigned int tdim(void) const 
      { 
         return tensorFeatureDimension; 
      }
      
      /** No bias
       */
      virtual bool bias() const
      {
         return false;
      }

      virtual std::vector<int>& Getcvmark()
      {
         return cvmark;
      }
};

#endif
