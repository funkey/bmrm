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

#ifndef _SEQFEATURE_CPP_
#define _SEQFEATURE_CPP_

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include "common.hpp"
#include "sml.hpp"
#include "configuration.hpp"
#include "bmrmexception.hpp"
#include "timer.hpp"
#include "seqfeature.hpp"

   
/**  Constructor
 */
CSeqFeature::CSeqFeature()   
        :seqfeature_verbosity(0),
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
         //density(0),
         X(0)       
{
        CTimer loadfeaturetime;

        // decide the format string to use
        if(sizeof(double) == sizeof(double)) 
                svec_feature_index_and_value_format = "%d:%lf\t";
        else 
                svec_feature_index_and_value_format = "%d:%f\t";
   
        // get configurations
        Configuration &config = Configuration::GetInstance();
   
        featureFile = config.GetString("Data.featureFile");
        if(config.IsSet("Data.verbosity"))
                seqfeature_verbosity = config.GetInt("Data.verbosity");   
   
        // read dataset into memory
        if(seqfeature_verbosity >= 1)
                std::cout << "Loading feature file... "<< std::endl;
   
        loadfeaturetime.Start();
        LoadFeatures();
        loadfeaturetime.Stop();
          
        // in centralized dataset, serial computation mode
        numOfAllSeq = numOfSeq;
        
        if(seqfeature_verbosity >= 2)
        {           
                std::cout << "loadfeaturetime : " << loadfeaturetime.CPUTotal() << std::endl;
        }        
}

CSeqFeature::~CSeqFeature()
{
        //if(X) delete[] X;
}


///** Read examples into memory
// */
//void CSeqFeature::LoadFeatures()
//{ 
//        unsigned int tmpFidx = 0;
//        double tmpFval = 0;
//        unsigned int featureCnt = 0;
//        unsigned int tmpSeqNum = 0, prevTmpSeqNum = 0;
//        bool firstSeq = true;
//        unsigned int phiNum = 0;
//        unsigned int posNum1 = 0, posNum2 = 0;
//        bool firstPhi2Pos1 = true, firstPhi2Pos2 = true;
//        std::string line = "";
//        std::string token = "";
//        std::ifstream featureFp;
//   
//        featureFp.open(featureFile.c_str());   
//        if(!featureFp.good()) 
//        {
//                string msg = "Cannot open feature file <" + featureFile + ">!";
//                throw CBMRMException(msg, "CSeqFeature::ScanFeatureFile()");
//        }
//   
//        // read header information
//        int headerInfoCnt = 3; // min duration, max duration, feature dimension
//        do {
//                getline(featureFp, line);
//                trim(line);
//                if(IsBlankLine(line)) continue;  // blank line
//                if(line[0] == '#') continue;  // comment line
//                if(sscanf(line.c_str(),"maxDuration:%d",&maxDuration)==1) headerInfoCnt--;
//                if(sscanf(line.c_str(),"minDuration:%d",&minDuration)==1) headerInfoCnt--;
//                if(sscanf(line.c_str(),"globalFeatureDim:%d",&featureDimension)==1) headerInfoCnt--;
//        } while(!featureFp.eof() && (headerInfoCnt != 0));
//        
//        assert(maxDuration >= minDuration);
//        assert(featureDimension < (1<<30));  // featureDimension is normally less then 1 billion
//                
//        if(featureFp.eof())
//                throw CBMRMException("Feature file does not contain valid examples","CSeqFeature::LoadFeatures()");
//        
//        // read sequences
//        nnz = 0;
//        seqfeature_struct tmp_seq;
//        vector<TheMatrix> tmp_phi_1;
//        vector<vector<TheMatrix> > tmp_phi_2;
//        while(!featureFp.eof()) 
//        {
//                // read sequence number
//                do {
//                        getline(featureFp, line);
//                        trim(line);
//                        if(IsBlankLine(line)) continue;  // blank line
//                        if(line[0] == '#') continue;  // comment line
//                        if(sscanf(line.c_str(),"sequence:%d",&tmpSeqNum)==1) break;
//                } while(!featureFp.eof());
//                
//                if(featureFp.eof())
//                        throw CBMRMException("Feature file does not contain valid phi:*","CSeqFeature::LoadFeatures()");
//                
//                if(!firstSeq && ((tmpSeqNum - prevTmpSeqNum) != 1))
//                        throw CBMRMException("sequence numbers must be consecutive and in increasing order","CSeqFeature::LoadFeatures()");
//                
//                prevTmpSeqNum = tmpSeqNum;
//                
//                // read phi:1 tag
//                phiNum = 0;
//                do {
//                        getline(featureFp, line);
//                        trim(line);
//                        if(IsBlankLine(line)) continue;  // blank line
//                        if(line[0] == '#') continue;  // comment line
//                        if(sscanf(line.c_str(),"phi:%d",&phiNum)==1) break;
//                } while(!featureFp.eof());
//                
//                if(featureFp.eof() || (phiNum != 1))
//                        throw CBMRMException("Feature file does not contain valid phi:1 tag","CSeqFeature::LoadFeatures()");                        
//                
//                // read phi:1 sparse vectors
//                do {
//                        getline(featureFp, line);
//                        trim(line);
//                        if(IsBlankLine(line)) continue;  // blank line
//                        if(line[0] == '#') continue;  // comment line
//                                                
//                        if(sscanf(line.c_str(),"phi:%d",&phiNum)==1)
//                                break;
//                        
//                        istringstream iss(line);
//                        iss >> token;
//                        if((sscanf(token.c_str(),"pos:%d",&posNum1) != 1))
//                        {
//                                cout << token << endl;
//                                throw CBMRMException("Feature file does not contain valid pos tag in phi:1","CSeqFeature::LoadFeatures()");
//                        }
//                        
//                        TheMatrix svec(1,featureDimension,SML::SPARSE);
//                        featureCnt = 0;
//                        while(!iss.eof())
//                        {
//                                iss >> token;
//                                if(sscanf(token.c_str(),svec_feature_index_and_value_format.c_str(),&tmpFidx, &tmpFval) != 2)
//                                {
//                                        ostringstream msg;
//                                        msg << "Invalid #" << featureCnt + 1 << " sparse vector element in phi:"<< phiNum << " seq:" << tmpSeqNum << " pos:" << posNum1;
//                                        throw CBMRMException(msg.str(),"CSeqFeature::LoadFeatures()");
//                                }
//                                svec.Set(0,tmpFidx,tmpFval);       
//                                nnz++;
//                                featureCnt++;
//                        }                                                
//                        phi_1.push_back(svec);
//                        
//                } while(!featureFp.eof());
//                
//                if(phi_1.size() < 1)
//                        throw CBMRMException("Feature file does not contain valid phi:1","CSeqFeature::LoadFeatures()");
//                
//                if(featureFp.eof() || (phiNum != 2))
//                        throw CBMRMException("Feature file does not contain valid phi:2","CSeqFeature::LoadFeatures()");
//                
//                // read phi:2 sparse vectors
//                unsigned int prevPosNum1 = 0, prevPosNum2 = 0; 
//                vector<TheMatrix> tmp_phi_2_svecs;
//                featureCnt = 0;
//                do {
//                        getline(featureFp, line);
//                        trim(line);
//                        if(IsBlankLine(line)) continue;  // blank line
//                        if(line[0] == '#') continue;  // comment line
//                        
//                        if((sscanf(line.c_str(),"phi:%d",&phiNum) == 1))
//                                break;
//                        
//                        istringstream iss(line);
//                        iss >> token;
//                        if((sscanf(token.c_str(),"pos:%d,%d",&posNum1,&posNum2) != 2))
//                                throw CBMRMException("Feature file contains invalid pos tag in phi:2","CSeqFeature::LoadFeatures()");
//                        
//                        if(prevPosNum2 != 0 && posNum2 !=0)
//                                firstPhi2Pos2 = false;
//                                
//                        if(prevPosNum1 !=0 && posNum1 != 0)
//                                firstPhi2Pos1 = false;
//                        
//                        if(!firstPhi2Pos2 && (prevPosNum1 == posNum1) && (prevPosNum2 >= posNum2))
//                        {
//                                ostringstream msg;
//                                msg << "previous posNum2 (" << prevPosNum2 << ") must be < current posNum2 (" << posNum2 << ") in phi:2 (phi:2 pos:" << posNum1 << "," << posNum2 << ")";
//                                throw CBMRMException(msg.str(),"CSeqFeature::LoadFeatures()");
//                        }
//
//                        
//                        if(!firstPhi2Pos1 && (prevPosNum1 >= posNum1))
//                        {
//                                ostringstream msg;
//                                msg << "previous posNum1 must be > current posNum1 in phi:2 (phi:2 pos:" << posNum1 << "," << posNum2 << ")";
//                                throw CBMRMException(msg.str(),"CSeqFeature::LoadFeatures()");
//                        }
//                        
//                        if(!firstPhi2Pos1 && (posNum1 != prevPosNum1))
//                        {
//                                phi_2.push_back(tmp_phi_2_svecs);
//                                tmp_phi_2_svecs.clear();                                
//                        }
//                        
//                        TheMatrix svec(1,featureDimension,SML::SPARSE);
//                        featureCnt = 0;
//                        while(!iss.eof())
//                        {
//                                iss >> token;
//                                if(sscanf(token.c_str(),svec_feature_index_and_value_format.c_str(),&tmpFidx, &tmpFval) != 2)
//                                {
//                                        ostringstream msg;
//                                        msg << "Invalid #" << featureCnt + 1 << " sparse vector element in seq:"<< tmpSeqNum << " phi:" << phiNum << " pos:" << posNum1;
//                                        throw CBMRMException(msg.str(),"CSeqFeature::LoadFeatures()");
//                                }
//                                svec.Set(0,tmpFidx,tmpFval);    
//                                nnz++;
//                                featureCnt++;
//                        }
//                        
//                        if(featureCnt == 0)
//                                throw CBMRMException("Feature file contains invalid phi:2 sparse vector","CSeqFeature::LoadFeatures()");
//                        
//                        tmp_phi_2_svecs.push_back(svec);
//                        prevPosNum2 = posNum2;
//                        prevPosNum1 = posNum1;
//
//                } while(!featureFp.eof());
//                
//                if(tmp_phi_2_svecs.size() > 0)
//                        phi_2.push_back(tmp_phi_2_svecs);
//                
//                if(phi_2.size() < 1)
//                        throw CBMRMException("Feature file contains no phi:2","CSeqFeature::LoadFeatures()");
//        }
//        
//        // data matrix density
//        density = ((double)nnz/featureDimension)/numOfSeq;
//   
//        featureFp.close();
//}


/** Read features into memory.
 *  
 *  Note:
 *  1. No serious (format) error checking.
 *  2. No comment lines allowed
 *
 */
void CSeqFeature::LoadFeatures()
{
        using namespace std;
        //printf("----------fuck------------\n");fflush(stdout);
        unsigned int tmpSeqID = 0;
        unsigned int tmpSeqSz = 0;
        unsigned int tmpPosNum1 = 0;
        unsigned int tmpPosNum2 = 0;
        unsigned int fidx = 0;
        double fval = 0.0;
        string line = "";
        string token = "";
        
//        FILE* ifp = fopen(featureFile.c_str(),"r");
//        fscanf(ifp,"#This dataset has %d sequences.\n",&numOfSeq);
//	fscanf(ifp,"maxDuration:%d\n",&maxDuration);
//	fscanf(ifp,"minDuration:%d\n",&minDuration);
//	fscanf(ifp,"globalFeatureDim:%d\n",&featureDimension);q
        
        ifstream ifp(featureFile.c_str());
        
        if(!ifp.good())
        {
                string msg = "Cannot open feature file <" + featureFile + ">!";
                throw CBMRMException(msg, "CSeqFeature::LoadFeatures()");  
        }
        
        // read header information
                // read header information
        int headerInfoCnt = 4; // num of sequence, min duration, max duration, feature dimension
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
        } while(!ifp.eof() && (headerInfoCnt != 0));

        if(seqfeature_verbosity >= 2)
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
                throw CBMRMException("Feature file does not contain valid examples","CSeqFeature::LoadFeatures()");
        
        //std::cout << "resizing X" << std::endl;
        //std::cout.flush();
        X.resize(numOfSeq);
        //std::cout << " done" << std::endl;
        //std::cout.flush();
        
        //std::cout << "creating tmp vector of thematrix for phi2..." << std::endl;
        //std::cout.flush();
        
        // pre-allocate vector of TheMatrix for use in creating phi2
        vector<TheMatrix*> tmp;
        //tmp.resize((maxDuration-minDuration)+1, TheMatrix(1,featureDimension,SML::SPARSE));
        
        //std::cout << "done" << std::endl;
        //std::cout.flush();
        
        for(unsigned int i=0;i < numOfSeq; i++)
        {                
                getline(ifp,line);
                trim(line);
                assert(! IsBlankLine(line));
                sscanf(line.c_str(),"sequenceID:%u\tsequenceSize:%u\n", &tmpSeqID, &tmpSeqSz);
                //fscanf(ifp,"sequenceID:%d\tsequenceSize:%d\n", &tmpSeqID, &tmpSeqSz);
                X[i].ID = tmpSeqID;
                X[i].len = tmpSeqSz;
                maxSeqLen = max(maxSeqLen, tmpSeqSz);
                minSeqLen = min(minSeqLen, tmpSeqSz);
                totalSeqLen += tmpSeqSz;
                
                if(seqfeature_verbosity >= 2)
                {
                        cout << "sequenceID:" << tmpSeqID << "   sequenceSize:"<< tmpSeqSz << endl;
                }
                                        
                // allocate space for phi fvec
                //std::cout << "resizing x phi1" << std::endl;
                //std::cout.flush();
                //X[i].phi_1.resize(tmpSeqSz, TheMatrix(1,featureDimension,SML::SPARSE));
                //std::cout << "done" << std::endl;
                //std::cout.flush();
                for(unsigned int j=0; j < tmpSeqSz; j++)
                {
                        //std::cout << "adding a vector... "; 
                        //std::cout.flush();
                        X[i].phi_1.push_back(new TheMatrix(1,featureDimension,SML::SPARSE));
                        //std::cout << "  done." << std::endl;
                        //std::cout.flush();
                }
                
                // allocate space for phi_2 fvec
                tmp.clear();
                for(unsigned int k=0; k<(maxDuration-minDuration)+1; k++)
                        tmp.push_back(new TheMatrix(1,featureDimension,SML::SPARSE));
                
                for(unsigned int j=0; j < tmpSeqSz; j++)                
                {
                        //std::cout << "adding a vector of vector ... "; 
                        //std::cout.flush();
                        X[i].phi_2.push_back(tmp);
                        //std::cout << "  done." << std::endl;
                        //std::cout.flush();
                }
                 
                // read phi:1 tag
                getline(ifp,line);
                trim(line);
                assert(!IsBlankLine(line));
                //fscanf(ifp,"phi:1\n");
                
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
                        //fscanf(ifp,"pos:%d\t",&tmpPosNum1);
                        //printf("pos:%d  ",tmpPosNum1); fflush(stdout);
                        //while(fscanf(ifp,svec_feature_index_and_value_format.c_str(),&fidx, &fval) == 2)
                        while(!iss.eof())
                        {               
                                iss >> token;
                                sscanf(token.c_str(),svec_feature_index_and_value_format.c_str(),&fidx, &fval);
                                X[i].phi_1[tmpPosNum1]->Set(0, fidx, fval);
                                nnz++;
                                //printf(" %d:%lf ",fidx,fval); fflush(stdout);
                        }
                        
                }
                
                // read phi:2 tag
                getline(ifp,line);
                trim(line);
                assert(!IsBlankLine(line));
                //fscanf(ifp,"phi:2\n");
                
                // read phi:2 svecs
                for(unsigned int j=0; j<tmpSeqSz; j++)
                {
                        //for(unsigned int k=minDuration; k < maxDuration && (j+k) <= tmpSeqSz; k++)
						for(unsigned int k=minDuration; k < maxDuration && (j+k) < tmpSeqSz; k++)
                        {
                                getline(ifp,line);
                                trim(line);
                                assert(!IsBlankLine(line));
                                istringstream iss(line);
                                iss >> token;
                                sscanf(token.c_str(),"pos:%u,%u", &tmpPosNum1, &tmpPosNum2);                                
                                //fscanf(ifp,"pos:%d,%d\t", &tmpPosNum1, &tmpPosNum2);                                
                                //printf("pos:%d,%d  ",tmpPosNum1,tmpPosNum2); fflush(stdout);
                                
                                //while(fscanf(ifp,svec_feature_index_and_value_format.c_str(),&fidx,&fval) == 2)
                                while(!iss.eof())
                                {
                                        iss >> token;
                                        sscanf(token.c_str(), svec_feature_index_and_value_format.c_str(), &fidx, &fval);
                                        X[i].phi_2[tmpPosNum1][tmpPosNum2-minDuration]->Set(0,fidx,fval);
                                        nnz++;
                                        //printf(" %d:%lf ",fidx,fval); fflush(stdout);
                                }
                        }
               }
       }

//       fclose(ifp);
        ifp.close();
}


void CSeqFeature::Dump()
{
        using namespace std;
        
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
                        for(unsigned int k=minDuration; k < maxDuration && (j+k) <= X[i].len; k++)
                        {
                                fprintf(ofp,"pos:%d,%d\t",j,k);
                                unsigned int len = X[i].phi_2[j][k-minDuration]->Length();
                                for(unsigned int n=0; n<len; n++)
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
