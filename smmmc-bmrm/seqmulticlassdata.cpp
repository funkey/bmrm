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
 * Last Updated: (28/10/2008))
 */

#ifndef _SEQMULTICLASSDATA_CPP_
#define _SEQMULTICLASSDATA_CPP_

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
#include "seqmulticlassdata.hpp"

using namespace std;

/**  Constructor
 */
CSeqMulticlassData::CSeqMulticlassData()
   : CData(), 
     CSeqMulticlassFeature(), 
     CSeqMulticlassLabel()
{
   // sanity check
   if(numOfAllSeq != numOfAllLabel)
   {
      ostringstream msg;
      msg << "Number of sequences (" << numOfAllSeq << ") does not match number of labels (" << numOfAllLabel << ")";
      throw CBMRMException(msg.str(),"CSeqMulticlassData::CSeqMulticlassData()");
   }
   
   // dataset statistics
   if(verbosity)
   {
      cout << "Dataset properties:"  << endl;
      cout << "1.  Feature file              : " << featureFile << endl;
      cout << "2.  Label file                : " << labelFile << endl;
      cout << "3.  Number of sequences       : " << numOfAllSeq << endl;    
      cout << "4.  Number of nonzero features: " << nnz << endl;
      cout << "5.  Max sequence length       : " << maxSeqLen << endl;
      cout << "6.  Min sequence length       : " << minSeqLen << endl;
      cout << "7.  Total sequence length     : " << totalSeqLen << endl;
      cout << "8.  Avg. sequence length      : " << (double)totalSeqLen/numOfAllSeq << endl;
      cout << "9.  Dim. of feature vector    : " << featureDimension << endl; 
      cout << "10. Start index of Phi1       : " << startIndexPhi1 << endl; 
      cout << "11. Start index of Phi2       : " << startIndexPhi2 << endl; 
      cout << "12. Start index of Phi3       : " << startIndexPhi3 << endl; 
      cout << endl;
   }	
}

void CSeqMulticlassData::SetCrossValidationData(unsigned int numOfFold, unsigned int foldIndex)
{
   // set a vector to indicate the examples' type as follows:
   // 0 -- train, 1 -- test, 2 -- validation
   unsigned int n_test, n_valid, n_train;
   unsigned int stride;
   unsigned int i;
   vector<int> cvmarktmp (numOfAllSeq,SMM::TRAIN_DATA);
   cvmark = cvmarktmp;
   
   n_test = (unsigned int)ceil(double(numOfAllSeq)/(numOfFold*2));
   n_valid = numOfAllSeq/numOfFold - n_test;
   n_train = numOfAllSeq-(n_test+n_valid);
   stride = numOfAllSeq/numOfFold;
   
   //set testing data
   for(i=0;i<n_test;i++)
   {
      cvmark[i+foldIndex*stride] = SMM::TEST_DATA;
   }
    
    
   //set validation data
   int validStart = n_test;
   for(i=0; i < n_valid; i++)
   {
      cvmark[i+(validStart)+foldIndex*stride] = SMM::VALID_DATA;
   }
   
   if(verbosity>=2)
   {	
      cout <<"cvmark:"<<endl;
      for(i=0;i<numOfAllSeq;i++)
         cout <<"["<<i<<"]:"<<cvmark[i]<<endl;
   }
}

/** Here tensorFeatureDimension is defined
 */
void CSeqMulticlassData::ExtendFeatures()
{
   unsigned int stepSize1 = startIndexPhi2 - startIndexPhi1;
   unsigned int stepSize2 = featureDimension-startIndexPhi2;   
   unsigned int numOfClass = this->getNumOfClass();
   unsigned int numOfPerson = this->getNumOfPerson();
   //extend to multi-class
   featureDimensionMid = 2+(stepSize1+1)*(numOfClass)+(stepSize2)*(numOfClass)+(stepSize2)*(numOfClass)*numOfClass;
   
   //extend to multi-person
   tensorFeatureDimension = featureDimensionMid*numOfPerson;
   //printf("fdim:%d pnum:%d cnum:%d tdim:%d\n ",featureDimensionMid,numOfPerson,numOfClass,tensorFeatureDimension);
}


void CSeqMulticlassData::TensorPhi1(TheMatrix* phi1, unsigned int classID, unsigned int personID, TheMatrix *v)
{
   v->Zero();
   assert(numOfClass > classID);
   
   unsigned int stepSize1 = startIndexPhi2-startIndexPhi1;	
   unsigned int lowBound = startIndexPhi1+(stepSize1+1)*(classID)+featureDimensionMid*personID;
   unsigned int upBound = startIndexPhi1+(stepSize1+1)*(classID+1)+featureDimensionMid*personID;
   unsigned int indexNew;
   
   //ToDO: should iterate non-zero sparse element, not iterate all dimensions
   //But how to do this in TheMatrix
   for(unsigned int m=0; m < featureDimension; m++)
   {	
      indexNew = m +(stepSize1+1)*(classID)+1+featureDimensionMid*personID;
      if((indexNew >= lowBound) && (indexNew < upBound))
      {			
         double val;
         phi1->Get(0,m,val);
         v->Set(0,indexNew,val);
      }
   }
   v->Set(0,1,1);
}

void CSeqMulticlassData::TensorPhi2(const TheMatrix*  phi2, unsigned int classIDPrev, unsigned int classID, unsigned int personID, int verbosity, TheMatrix *v)
{
   v->Zero();    
   assert(numOfClass > classID);
   
   unsigned int stepSize1 = startIndexPhi2-startIndexPhi1;	
   unsigned int stepSize2 = featureDimension-startIndexPhi2;
   unsigned int indexNew2, indexNew3, lowBound2, lowBound3, upBound2, upBound3;
   
   //get dim index bounds
   //lowBound2 = startIndexPhi1 + (stepSize1+1)*(numOfClass) + stepSize2*(classIDPrev) + featureDimensionMid*personID;
   //upBound2  = startIndexPhi1 + (stepSize1+1)*(numOfClass) + stepSize2*(classIDPrev+1) + featureDimensionMid*personID;		
   //lowBound3 = startIndexPhi1 + (stepSize1+1)*(numOfClass) + stepSize2*(numOfClass) + stepSize2*(numOfClass*classIDPrev+classID) + featureDimensionMid*personID;
   //upBound3  = startIndexPhi1 + (stepSize1+1)*(numOfClass) + stepSize2*(numOfClass) + stepSize2*(numOfClass*classIDPrev+classID+1) + featureDimensionMid*personID;
   
   lowBound2 = startIndexPhi1 + (stepSize1+1)*(numOfClass) + featureDimensionMid*personID + stepSize2*(classIDPrev);
   upBound2  = lowBound2 + stepSize2;
   lowBound3 = startIndexPhi1 + (stepSize1+1)*(numOfClass) + stepSize2*(numOfClass) + featureDimensionMid*personID + stepSize2*(numOfClass*classIDPrev+classID);
   upBound3  = lowBound3 + stepSize2;
   
   if(verbosity >= 1)
   {
      printf("low2:%d up2:%d low3:%d up3:%d",lowBound2,upBound2,lowBound3,upBound3);
      phi2->Print();
      fflush(stdout);
   }
   
   for(unsigned int m=0; m < featureDimension; m++)
   {       
      indexNew2 = 0;
      indexNew3 = 0;		
      
      if(m < startIndexPhi3)
      {
         indexNew2 = (m-startIndexPhi2) + 2 + (stepSize1+1)*(numOfClass) + stepSize2*(classIDPrev) + featureDimensionMid*personID;
      }
      else 
         indexNew3 = (m-startIndexPhi2) + 2 + (stepSize1+1)*(numOfClass) + stepSize2*(numOfClass) + stepSize2*(numOfClass*classIDPrev+classID) + featureDimensionMid*personID;
      
      // extend phi2 to phi2
      if((indexNew2 >= lowBound2) && (indexNew2 < upBound2))
      {
         double val;
         phi2->Get(0,m,val);			
         v->Set(0,indexNew2,val);
      }
      
      // extend phi2 to phi3
      if((indexNew3 >= lowBound3) && (indexNew3 < upBound3))
      {
         double val;
         phi2->Get(0,m,val);
         v->Set(0,indexNew3,val);
      }
   }
}

#endif
