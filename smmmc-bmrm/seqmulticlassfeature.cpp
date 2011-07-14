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

#ifndef _SEQMULTICLASSFEATURE_CPP_
#define _SEQMULTICLASSFEATURE_CPP_

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include "common.hpp"
#include "sml.hpp"
#include "configuration.hpp"
#include "bmrmexception.hpp"
#include "seqmulticlassfeature.hpp"

using namespace std;
   
/**  Constructor
 */
CSeqMulticlassFeature::CSeqMulticlassFeature()   
    :seqmulticlassfeature_verbosity(0),
     maxDuration(0),
     minDuration(1<<30),
     numOfSeq(0),
     numOfAllSeq(0),                
     featureDimension(0),                
     maxSeqLen(0),
     minSeqLen(1<<30),
     totalSeqLen(0),
     featureFile(""),
     nnz(0),
     X(0)       
{
    // get configurations
    Configuration &config = Configuration::GetInstance();
    
    featureFile = config.GetString("Data.featureFile");
    if(config.IsSet("Data.verbosity"))
	seqmulticlassfeature_verbosity = config.GetInt("Data.verbosity");   
    
    // read dataset into memory
    if(seqmulticlassfeature_verbosity >= 1)
	cout << "Loading feature file... "<< endl;
    
    LoadFeatures();
          
    // in centralized dataset, serial computation mode
    numOfAllSeq = numOfSeq;
}

CSeqMulticlassFeature::~CSeqMulticlassFeature()
{
    //if(X) delete[] X;
}


/** Read features into memory.
 *  
 *  Note:
 *  1. No serious (format) error checking.
 *  2. No comment lines allowed
 *
 */
void CSeqMulticlassFeature::LoadFeatures()
{
   unsigned int tmpSeqID = 0;
   unsigned int tmpSeqSz = 0;
   unsigned int tmpPosNum1 = 0;
   unsigned int tmpPosNum2 = 0;
   unsigned int fidx = 0;
   double fval = 0.0;
   string line = "";
   string token = "";
   
   ifstream ifp(featureFile.c_str());
   
   if(!ifp.good())
   {
      string msg = "Cannot open feature file <" + featureFile + ">!";
      throw CBMRMException(msg, "CSeqMulticlassFeature::LoadFeatures()");  
   }
   
   
    // read header information
   int headerInfoCnt = 7; // num of sequence, min duration, max duration, feature dimension,startIndexPhi1,startIndexPhi2,startIndexPhi3
   do {
      getline(ifp, line);
      //cout << line.c_str() << endl;fflush(stdout);
      
      trim(line);
      
      if(IsBlankLine(line)) continue;  // blank line
      if(line[0] == '#') continue;  // comment line
      if(sscanf(line.c_str(),"numOfSequences:%u",&numOfSeq)==1) headerInfoCnt--;
      if(sscanf(line.c_str(),"maxDuration:%u",&maxDuration)==1) headerInfoCnt--;
      if(sscanf(line.c_str(),"minDuration:%u",&minDuration)==1) headerInfoCnt--;
      if(sscanf(line.c_str(),"globalFeatureDim:%u",&featureDimension)==1) headerInfoCnt--;
      if(sscanf(line.c_str(),"startIndexPhi1:%u",&startIndexPhi1)==1) headerInfoCnt--;
      if(sscanf(line.c_str(),"startIndexPhi2:%u",&startIndexPhi2)==1) headerInfoCnt--;
      if(sscanf(line.c_str(),"startIndexPhi3:%u",&startIndexPhi3)==1) headerInfoCnt--;
   } while(!ifp.eof() && (headerInfoCnt != 0));
   
   if(seqmulticlassfeature_verbosity >= 2)
   {
      cout << "No. of sequences: " << numOfSeq << endl;
      cout << "Max duration    : " << maxDuration << endl;
      cout << "Min duration    : " << minDuration << endl;                
      cout << "Feature dim.    : " << featureDimension << endl;
   }
   
   
   assert(numOfSeq < (1<<30)); 
   assert(maxDuration >= minDuration);
   assert(featureDimension < (1<<30));  // featureDimension is normally less then 1 billion
   
   if(ifp.eof())
      throw CBMRMException("Feature file does not contain valid examples","CSeqMulticlassFeature::LoadFeatures()");
   
   X.resize(numOfSeq);
   
   for(unsigned int i=0;i < numOfSeq; i++)
   {
      if(i==0)
      {
         getline(ifp,line);
         trim(line);
         assert(! IsBlankLine(line));
         sscanf(line.c_str(),"sequenceID:%u\n", &tmpSeqID);
      }
      else
         tmpSeqID = i;
      
      getline(ifp,line);
      trim(line);
      assert(! IsBlankLine(line));
      sscanf(line.c_str(),"sequenceSize:%u\n", &tmpSeqSz);
      
      X[i].ID = tmpSeqID;
      X[i].len = tmpSeqSz;
      maxSeqLen = max(maxSeqLen, tmpSeqSz);
      minSeqLen = min(minSeqLen, tmpSeqSz);
      totalSeqLen += tmpSeqSz;
      
      if(seqmulticlassfeature_verbosity >= 2)
      {
         cout << "sequenceID:" << tmpSeqID << "   sequenceSize:"<< tmpSeqSz << endl;
      }
      
      // allocate space for phi fvec
      for(unsigned int j=0; j < tmpSeqSz; j++)
      {
         X[i].phi_1.push_back(new TheMatrix(1,featureDimension,SML::SPARSE));
      }
      
      // allocate space for phi_2 fvec
      for(unsigned int j=0; j < tmpSeqSz; j++)                
      {
         vector<TheMatrix*> tmp;
         tmp.clear();
         for(unsigned int k=0; k<(maxDuration-minDuration)+1; k++)
            tmp.push_back(new TheMatrix(1,featureDimension,SML::SPARSE));
         
         X[i].phi_2.push_back(tmp);
      }
      
      // read phi:1 tag
      getline(ifp,line);
      trim(line);
      assert(!IsBlankLine(line));
      
      // read phi:1 svecs
      for(unsigned int j=0; j<tmpSeqSz; j++)
      {
         TheMatrix sv(1,featureDimension, SML::SPARSE);
         getline(ifp,line);
         trim(line);
         assert(!IsBlankLine(line));
         istringstream iss(line);
         iss >> token;
         sscanf(token.c_str(),"pos:%u",&tmpPosNum1);
         
         X[i].phi_1[tmpPosNum1]->Zero();
         while(!iss.eof())
         {               
            iss >> token;
            sscanf(token.c_str(),"%u:%lf\t",&fidx, &fval);
            X[i].phi_1[tmpPosNum1]->Set(0, fidx, fval);
            nnz++;
         }
      }
      
      // read phi:2 tag
      getline(ifp,line);
      trim(line);
      assert(!IsBlankLine(line));
      
      // read phi:2 svecs
      while(! ifp.eof())
      {
         getline(ifp,line);
         if(line[0]!='p')
         {
            break;
         }
         trim(line);
         
         
         assert(!IsBlankLine(line));
         istringstream iss(line);
         iss >> token;
         sscanf(token.c_str(),"pos:%u,%u", &tmpPosNum1, &tmpPosNum2);   
         X[i].phi_2[tmpPosNum1][tmpPosNum2-minDuration]->Zero();
         while(!iss.eof())
         {
            iss >> token;
            sscanf(token.c_str(), "%u:%lf\t", &fidx, &fval);
            X[i].phi_2[tmpPosNum1][tmpPosNum2-minDuration]->Set(0,fidx,fval);
            nnz++;
         }
      }
   }
   ifp.close();
}


void CSeqMulticlassFeature::Dump()
{
   FILE *ofp = fopen("seqfeature.dump","w");
   fprintf(ofp,"numSequences:%u\n",numOfSeq);
   fprintf(ofp,"maxDuration:%u\n",maxDuration);
   fprintf(ofp,"minDuration:%u\n",minDuration);
   fprintf(ofp,"globalFeatureDim:%d\n",featureDimension);        
    
   for(unsigned int i=0; i<numOfSeq; i++)
   {
      fprintf(ofp,"sequenceID:%d\tsequenceSize:%d\n",X[i].ID, X[i].len);
      fprintf(ofp,"phi:1\n");
      for(unsigned int j=0; j< X[i].len; j++)
      {
         fprintf(ofp,"pos:%d\t",j);
         unsigned int len = X[i].phi_1[j]->Length();
	    for(unsigned int k=0; k<len; k++)
	    {
           double val = 0.0;
           X[i].phi_1[j]->Get(0,k,val);
           if(fabs(val) > 1e-20)
           {
              fprintf(ofp,"%d:%f\t",k,val);
           }
	    }
	    fprintf(ofp,"\n");
      }
      
      fprintf(ofp,"phi:2\n");
      for(unsigned int j=0; j < X[i].len; j++)
      {
         //for(unsigned int k=0; k < X[i].phi_2[j].size(); k++)
         //for(unsigned int k=0; k < (maxDuration-minDuration && ); k++)
         for(unsigned int k=minDuration; k < (maxDuration) && (j+k) <= X[i].len; k++)
         {
            fprintf(ofp,"pos:%d,%d\t",j,k);
            //printf("pos:%d,%d\t",j,k);
            //unsigned int len = X[i].phi_2[j][k-minDuration]->Length();
            //printSparseVector(featureDimension,X[i].phi_2[j][k-minDuration]);
            for(unsigned int n=0; n<featureDimension; n++)
            {
               double val = 0.0;
               X[i].phi_2[j][k-minDuration]->Get(0,n,val);
               if(fabs(val) > 1e-20)
               {
                  fprintf(ofp, "%d:%f\t",n,val);
               }
            }
            fprintf(ofp,"\n");
         }
      }
   }
   fclose(ofp);
}


#endif
