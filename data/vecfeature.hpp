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
 * Last Updated: 19/04/2008
 */

#ifndef _VECFEATURE_HPP_
#define _VECFEATURE_HPP_

#include <vector>
#include <iostream>

#include "common.hpp"
#include "sml.hpp"


/** Container for feature vectors
 *    
 *   Feature vector file format: rows of <feature_vector>
 *   where
 *   <feature_vector> .=. [qid:<id>] <dense_vector> | <sparse_vector>
 *   <dense_vector>   .=. <value_1> [<value_2> ... <value_k>]
 *   <sparse_vector>  .=. <index_1>:<value_1> [<index_2>:<value_2> ... <index_k>:<value_k>]
 *   <id>             .=. a number (examples with same qid must be adjacent to each other in the data file)
 *   <index_j>        .=. index (positive integer) of j-th non-zero element of a particular sparse feature vector
 *                        Note that index_i > index_j for j>i
 *   <value_j>        .=. value (scalar) of j-th non-zero element of a particular feature vector
 *
 * Notes:
 *   1.  The data loading starts during the object construction.
 *   2.  For centralized dataset and serial computation mode only!
 */
class CVecFeature
{
   protected:               
      /** verbosity level
       */
      int vecfeature_verbosity;
      
      /** Feature matrix (each row is an example)
       */
      TheMatrix* X;
      
      /** A list of examples as (explicit) row vectors
       */
      TheMatrix** x;
      
      /** Types of feature vector
       */
      enum FEATURE_TYPE {SPARSE_FEATURE, DENSE_FEATURE};
      
      /** whether feature vector comes with artificial bias feature
       */
      bool biasFlag;
      
      /** Index of the first example to load
       */
      unsigned int startExampleIndex;
      
      /** number of examples in this sub-dataset
       */
      unsigned int numOfExample;
      
      /** number of examples in the whole dataset
       */
      unsigned int numOfAllExample;
      
      /** dimensionality of the feature vector
       */
      unsigned int featureDimension;        
      
      /** Index of the first subset to load
       */
      unsigned int startSubsetIndex;
      
      /** Number of subsets in this (sub-)dataset
       */
      unsigned int numOfSubset; 
      
      /** Number of subsets in the dataset
       */
      unsigned int numOfAllSubset;
      
      /** Number of examples in each subset
       */
      std::vector<int>* subsetSizes; 
      
      /** To create matrix row view for feature matrix
       */
      bool featureMatrixRowView;
      
      /** Flag for type of feature vector e.g., dense or sparse
       */
      unsigned int featureType;
      
      /** Name of the file containing feature vectors
       */
      std::string featureFile;
      
      /** A string template for sparse vector element
       */
      std::string svec_feature_index_and_value_format;
      
      /** A string template for dense vector element
       */
      std::string scalar_value_format;
      
      /** the set of nnz per example
       */
      std::vector<int> nnzSet;
      
      /** number of nonzero features in whole dataset
       */
      unsigned int allnnz;
      
      /** maximum nnz per feature vector
       */
      unsigned int maxNNZ;
      
      /** number of nonzero features in sub-dataset
       */
      unsigned int nnz;
      
      /** nnz / total number of entries in feature matrix
       */
      double density;
      
      /** Structure to keep some info. of subsets of examples in the dataset
       */
      struct Subset {
            int ID;
            int startIndex;
            int size;
      };
      
      /** Scan the feature file and determine some feature set properties such as dimension of example        
       */
      virtual void ScanFeatureFile();
      
      /** Allocate data matrix and load fetures from data file.
       *  Read "numOfExample" starting from the "startExampleIndex".
       */
      virtual void LoadFeatures();     
      
   public:       
      /** Constructor
       */
      CVecFeature();
      
      
      /** Destructor
       */
      virtual ~CVecFeature();
      
      /** Load feature matrix.
       *  The total number of units N in the dataset will be divided 
       *  into "nparts" parts according to the following rule:
       *  size(part_i) = floor(N/nparts) + [i <= N%nparts], for i=1,...,nparts, where [a] is 1 if a is true, 0 otherwise
       *  Then load only i-th part into memory (0-based indexing)
       *  
       *  Note:
       *  1. when number of subsets > 1, 1 unit := 1 subset (a group of data points) 
       *     otherwise, 1 unit := 1 data point
       *  2. when start = nparts = 1, we load whole dataset
       
       */
      virtual void LoadFeatureMatrix(unsigned int start=0, unsigned int nparts=1);
      
      
      /** Given a weight vector w return the prediction f = Xw. It is the
       *  duty of the caller to ensure that w and f have conforming
       *  dimensions.  
       * 
       *  @param w [read]  the weight vector
       *  @param f [write] the prediction 
       */
      virtual void ComputeF(const TheMatrix& w, TheMatrix& f);
      
      
      /** Given a weight vector w return the prediction f_i = x_i*w. 
       *  It is the duty of the caller to ensure that w and f have conforming
       *  dimensions.  
       * 
       *  @param w [read]  the weight vector (i.e. a matrix)
       *  @param f [write] the prediction 
       */
      virtual void ComputeFi(const TheMatrix& w, TheMatrix& f, const unsigned int i);
      
      /** Given a weight vector w return the prediction f = x_i*w. 
       *  It is the duty of the caller to ensure that w has conforming
       *  dimensions.  
       * 
       *  @param w [read]  the weight vector
       */
      virtual void ComputeFi(const TheMatrix& w, double & f, const unsigned int i);
      
      
      /** f = X*w, where X is a data matrix (with rows of feature vectors) 
       */
      virtual void XMultW(const TheMatrix& w, TheMatrix& f)
      {
         X->Dot(w, f);
      }

      
      /** f = X^T*w, where X is a data matrix (with rows of feature vectors) 
       */
      virtual void XTMultW(const TheMatrix& w, TheMatrix& f)
      {
         X->TransposeDot(w, f);
      }
      
      
      /** f = w*X, where X is a data matrix (with rows of feature vectors) 
       */
      virtual void WMultX(const TheMatrix& w, TheMatrix& f)
      {
         w.Dot(*X, f);
      }
      
      
      /** f = w*X, where X is a data matrix (with rows of feature vectors) 
       */
      virtual void WTMultX(const TheMatrix& w, TheMatrix& f)
      {
         w.TransposeDot(*X, f);
      }
      
      /** Adds scale times x_i to w. 
       *  It is the duty of the caller to ensure that w has conforming
       *  dimensions.
       *
       *  @param w [write] weight vector
       *  @param i [read] position of the x_i
       *  @param scale [read] scaling factor    
       */
      virtual void AddElement(TheMatrix& w, const unsigned int &i,double scale);

      
      /** Given the loss vector l return the gradient grad = lX. It is the
       *  duty of the caller to ensure that l and grad have conforming
       *  dimensions.
       * 
       *  @param l    [read]  the loss vector
       *  @param grad [write] the gradient
       */
      virtual void Grad(const TheMatrix& l, TheMatrix& grad);
      
      
      /** Create explicit row views of feature matrix
       */
      virtual void CreateFeatureMatrixRowViews();  // create matrix row view for every example vector
      
      /** the index of the first subset to read
       */
      unsigned int StartSubSetIndex() {return startSubsetIndex;}
      
      /** return the number of subsets of this sub-dataset
       */
      unsigned int NumOfSubset() {return numOfSubset;}
      
      /** return the number of subsets of the WHOLE dataset when used in distributed environment
       *  for centralized dataset, this is equivalent to NumOfSubset()
       */
      unsigned int NumOfAllSubset() {return numOfAllSubset;}
      
      
      unsigned int StartExampleIndex() {return startExampleIndex;}
      unsigned int NumOfExample() {return numOfExample;}
      unsigned int NumOfAllExample() {return numOfAllExample;}
      bool HasMatrixRowView() {return featureMatrixRowView;}
      bool HasBiasFlag() {return biasFlag;}
      bool HasData() {return (X != 0);}
      unsigned int FeatureDimension() {return featureDimension;}
      unsigned int NNZ() {return nnz;}
      double Density() {return density;}
      std::string FeatureFile() {return featureFile;}
      
      void ShowFeature() { X->Print(); }
      
      /** subset information
       */
      Subset *subset;              
};

#endif
