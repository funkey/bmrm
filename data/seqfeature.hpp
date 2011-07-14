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
 * Created: (29/01/2008) 
 *
 * Last Updated:
 */

#ifndef _SEQFEATURE_HPP_
#define _SEQFEATURE_HPP_

#include <vector>
#include <iostream>

#include "common.hpp"
#include "sml.hpp"


/** Container for sequence feature vectors
 *  Mainly for Automatic Paragraph Segmentation 
 *  
 *   Example format: 
 *  
 *   maxDuration:<P>
 *   minDuration:<N>
 *   globalFeatureDim:<P>
 *   sequence:<N> 
 *   phi:1 
 *   pos:<N> <sparse_vector>
 *   ...
 *   phi:2 
 *   pos:<N>,<N> <sparse_vector>
 *   ...
 *   sequence:<N>
 *   phi:1 
 *   pos:<N> <sparse_vector>
 *   ...
 *   phi:2 
 *   pos:<N>,<N> <sparse_vector>
 *   ...
 *
 *   where
 *   <P>              .=. positive integer
 *   <N>              .=. natural number
 *   <sparse_vector>  .=. <index_1>:<value_1> [<index_2>:<value_2> ... <index_k>:<value_k>]
 *   <index_j>        .=. index (positive integer) of j-th non-zero element of a particular sparse feature vector
 *                        Note that index_i > index_j for j>i
 *   <value_j>        .=. value (scalar) of j-th non-zero element of a particular feature vector
 *
 * Notes:
 *   1.  The data loading starts during the object construction.
 *   2.  For centralized dataset and serial computation mode only!
 */
class CSeqFeature
{
public:
        /** sequence structure
         */
        struct seqfeature_struct {
                unsigned int ID;
                unsigned int len;
                std::vector<TheMatrix*> phi_1;
                std::vector<std::vector<TheMatrix*> > phi_2;
        };
        
        
        CSeqFeature();
        virtual ~CSeqFeature();
      
        virtual unsigned int MaxDuration() {return maxDuration;}
        virtual unsigned int MinDuration() {return minDuration;}
        virtual unsigned int MaxSequenceLength() {return maxSeqLen;}
        virtual unsigned int MinSequenceLength() {return minSeqLen;}
        virtual const std::vector<seqfeature_struct> &features() { return X;}        
        void Dump();
        
protected:               
        /** verbosity level
         */
        int seqfeature_verbosity;                
        
        /** Maximum duration of a segment in a sequence
         */
        unsigned int maxDuration;
        
        /** Minimum duration of a segment in a sequence
         */
        unsigned int minDuration;
        
        /** number of examples in this sub-dataset
         */
        unsigned int numOfSeq;

        /** number of examples in the whole dataset
         */
        unsigned int numOfAllSeq;

        /** dimensionality of the feature vector
         */
        unsigned int featureDimension;        

        /** length of longest sequence
         */
        unsigned int maxSeqLen;
        
        /** length of shortest sequence
         */
        unsigned int minSeqLen;
        
	/** total length of all sequences
	 */
	unsigned int totalSeqLen;

        /** Name of the file containing feature vectors
         */
        std::string featureFile;
      
        /** A string template for sparse vector element
         */
        std::string svec_feature_index_and_value_format;
      
        /** number of nonzero features
         */
        unsigned int nnz;

        /** numOfNonzero / total number of entries in feature matrix
         */
        //double density;
        
        /** Sequences
         */
        std::vector<seqfeature_struct> X;

               
        /** Allocate data matrix and load features from data file
         */
        virtual void LoadFeatures();     
        

};

#endif
