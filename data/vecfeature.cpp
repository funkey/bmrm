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

#ifndef _VECFEATURE_CPP_
#define _VECFEATURE_CPP_

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include "common.hpp"
#include "sml.hpp"
#include "configuration.hpp"
#include "bmrmexception.hpp"
#include "timer.hpp"
#include "vecfeature.hpp"

using namespace std;

   
/**  Constructor
 */
CVecFeature::CVecFeature()
    :vecfeature_verbosity(0),
     X(0),
     x(0),
     biasFlag(false),
     startExampleIndex(0),
     numOfExample(0),
     numOfAllExample(0),                
     featureDimension(0), 
     startSubsetIndex(0),
     numOfSubset(0),
     numOfAllSubset(0),
     subsetSizes(new vector<int>),
     featureMatrixRowView(false),
     featureType(SPARSE_FEATURE),
     featureFile(""),
     allnnz(0),
     maxNNZ(0),
     nnz(0),
     density(0),
     subset(0)       
{}


/** destructor
 **/
CVecFeature::~CVecFeature()
{ 
    if(X) delete X;
    if(x) 
    {
        for(unsigned int i=0; i<numOfExample; i++) 
            delete x[i];
        free(x); 
    }  
    
    if(subset) delete [] subset;
}



/** Load feature matrix
 *  Load only the "start"-th part of total "nparts" parts
 */
void CVecFeature::LoadFeatureMatrix(unsigned int start, unsigned nparts)
{
    assert(start <= nparts); 
    
    CTimer scanfeaturetime;
    CTimer loadfeaturetime;
    
    svec_feature_index_and_value_format = "%d:%lf";
    scalar_value_format = "%lf";
    
   
    // get configurations
    Configuration &config = Configuration::GetInstance();
        
    featureFile = config.GetString("Data.featureFile");
    if(config.IsSet("Data.verbosity"))
        vecfeature_verbosity = config.GetInt("Data.verbosity");
    
    if(config.IsSet("Data.bias"))
        biasFlag = config.GetBool("Data.bias");
    
    // collect some properties of the dataset
    if(vecfeature_verbosity >= 1)
        cout << "Scanning feature file... " << endl;
    
    scanfeaturetime.Start();
    ScanFeatureFile();
    scanfeaturetime.Stop();
    
    
    // compute number of examples and subsets       
    // 1 unit := 1 subset
    if(numOfAllSubset > 1)
    {                
        unsigned int q = numOfAllSubset / nparts;
        unsigned int r = numOfAllSubset % nparts;
        numOfSubset = q + ((start < r) ? 1 : 0);        
        startSubsetIndex = start * q + ((start<r) ? start : r);
        
        //printf("startsubsetindex=%d\n",startSubsetIndex);
        
        // subset information
        unsigned int exIdx = 0;                
        for(unsigned int i=0; i<startSubsetIndex; i++)
            exIdx += subsetSizes->at(i);
        
        startExampleIndex = exIdx;                                
        subset = new Subset[numOfSubset];                
        numOfExample = 0;
        for(unsigned int i=0; i<numOfSubset; i++) 
        {
            subset[i].startIndex = exIdx;
            subset[i].size = subsetSizes->at(startSubsetIndex + i);
            subset[i].ID = 0;
            exIdx += subset[i].size;
            numOfExample += subset[i].size;                        
        }
    }
    // 1 unit := 1 point
    else
    {
        numOfSubset = 1;
        startSubsetIndex = 0;
        
        unsigned int q = numOfAllExample / nparts;
        unsigned int r = numOfAllExample % nparts;
        numOfExample = q + ((start < r) ? 1 : 0);        
        startExampleIndex = start * q + ((start<r) ? start : r);
        
        // subset information
        subset = new Subset[1];
        subset[0].startIndex = 0;
        subset[0].size = subsetSizes->at(0);
        subset[0].ID = 0;                
    }
        
    assert(numOfExample > 0);
    assert(numOfSubset > 0);
    assert(startExampleIndex >= 0);
    assert(numOfAllExample == nnzSet.size());
    
    // compute nnz for sub-dataset
    for(unsigned int i=0; i<numOfExample; i++)
        nnz += nnzSet[startExampleIndex + i];
    
    // compute feature matrix density
    density = ((double)nnz/featureDimension)/numOfExample;
    
    // read dataset into memory
    if(vecfeature_verbosity >= 1)
        cout << "Loading feature file... "<< endl;
    
    loadfeaturetime.Start();
    LoadFeatures();
    loadfeaturetime.Stop();
    
    if(vecfeature_verbosity >= 2)
    {           
        cout << "scanfeaturetime : " << scanfeaturetime.CPUTotal() << endl;
        cout << "loadfeaturetime : " << loadfeaturetime.CPUTotal() << endl;
    }                
}


/**  Determine the following properties of the feature vector file:
 *   1.  dimensionality of example
 *   2.  maximum of number of nonzero in feature vectors
 *   3.  number of examples
 *   4.  density of data matrix
 *   5.  number of nonzero feature elements
 *   6.  type of feature vector
 *   
 */
void CVecFeature::ScanFeatureFile()
{ 
    unsigned int prevFeatureDimension = 0;
    unsigned int prevFidx = 0;
    unsigned int tmpFidx = 0;
    double tmpFval = 0;
    unsigned int lineCnt = 0;
    unsigned int featureCnt = 0;
    bool firstQuery = true;
    unsigned int prevQID = 0;
    unsigned int curQID = 0;
    string line = "";
    string token = "";
    ifstream featureFp;
    
    featureFp.open(featureFile.c_str());   
    if(!featureFp.good()) 
    {
        string msg = "Cannot open feature file <" + featureFile + ">!";
        throw CBMRMException(msg, "CVecFeature::ScanFeatureFile()");
    }
    
    // Read first line to check whether data file is in SPARSE or DENSE format
    bool done = false;
    while(!done && !featureFp.eof()) 
    {
        getline(featureFp, line);    
        trim(line);
        if(IsBlankLine(line)) continue;
        if(line[0] == '#') continue;  // comment line
        
        istringstream iss(line);
        while(!iss.eof()) 
        {
            iss >> token;
            if(token[0] == '#') break;
            if(token.find("id:") != string::npos) continue;
            
            if(sscanf(token.c_str(), svec_feature_index_and_value_format.c_str(), &tmpFidx, &tmpFval) == 2)
            {
                featureType = SPARSE_FEATURE;
                if(vecfeature_verbosity)
                    cout << "Sparse feature file detected... " << endl;
                done = true;
                break;
            }
            else 
            {        
                featureType = DENSE_FEATURE;
                if(vecfeature_verbosity) 
                    cout << "Dense data file detected... " << endl;
                done = true;
                break;
            }
        }
    }
    
   
    // do a pre-scan base on the feature type detected
    line = "";
    featureFp.seekg(ios::beg);
    
    if(featureType == SPARSE_FEATURE)
    {
        for(numOfAllExample=0, lineCnt=0; !featureFp.eof(); lineCnt++) 
        {
            getline(featureFp, line);
            trim(line);
            if(IsBlankLine(line)) continue;  // blank line
            if(line[0] == '#') continue;  // comment line
            istringstream iss(line);
            
            // collect information
            // query id is always put in front of features !
            size_t pos = line.find("qid:");
            if(pos != string::npos) 
            {
                sscanf(&line[pos], "qid:%u ", &curQID);      
                iss >> token; // throw this token away
            }
            
            // keep track of number of subset (e.g., query)
            if(firstQuery || (prevQID != curQID))
            {
                firstQuery = false;
//                                // check if qid's are in ascending order
//                                if(prevQID > curQID) 
//                                {
//                                        ostringstream ostr;
//                                        ostr << "QID's of line #" << lineCnt-1 << " and line #" << lineCnt << " are out of order!\n";
//                                        throw CBMRMException(ostr.str(), "CVecFeature::ScanFeatureFile()");
//                                }
            
                numOfAllSubset++;
                subsetSizes->push_back(0);
                prevQID = curQID;
            }
            
            // feature pairs
            for(featureCnt=0; !iss.eof();)
            {    
                iss >> token;            
                if(IsBlankLine(token)) break;
                if(token[0] == '#') break;
                
                if(sscanf(token.c_str(), svec_feature_index_and_value_format.c_str(), &tmpFidx, &tmpFval) < 2) 
                {
                    ostringstream ostr;
                    ostr << "Feature file contains invalid feature pair at line " << lineCnt+1 << "!\n"
                         << "Token: " << token << "\n";
                    throw CBMRMException(ostr.str(), "CVecFeature::ScanFeatureFile()");
                }
                
                // check fidx order
                if(prevFidx >= tmpFidx) 
                {
                    ostringstream ostr;
                    ostr << "Feature file contains feature indices at line " << lineCnt+1 << "!\n"
                         << "Line: " << line << "\n";
                    throw CBMRMException(ostr.str(), "CVecFeature::ScanFeatureFile()");
                }
                
                // update statistics
                featureCnt++;
                featureDimension = max(featureDimension, tmpFidx);
            }
            
            if(featureCnt <= 0) 
            {
                ostringstream ostr;
                ostr << "Example #" << numOfAllExample+1 << " has no non-zero feature element!\n";
                throw CBMRMException(ostr.str(), "CVecFeature::ScanFeatureFile()");
            }
            
            // update statistics
            if(biasFlag)
                featureCnt += 1;    // for shifted hyperplane
            subsetSizes->at(numOfAllSubset-1) += 1;
            nnzSet.push_back(featureCnt);
            allnnz += featureCnt;
            maxNNZ = max(maxNNZ, featureCnt);
            numOfAllExample++;
        }
    }
    else if(featureType == DENSE_FEATURE)
    {     
        for(numOfAllExample=0, lineCnt=0; !featureFp.eof(); lineCnt++)
        {
            getline(featureFp, line);
            trim(line);
            if(IsBlankLine(line)) continue;  // blank line
            if(line[0] == '#') continue;  // comment line
            istringstream iss(line);
            
            // collect information
            // query id is always in front of feature vector!
            size_t pos = line.find("qid:");
            if(pos != string::npos) 
            {
                sscanf(&line[pos], "qid:%u ", &curQID);      
                iss >> token; // throw away this token
            }
            
            // keep track of number of subset (e.g., query)
            if(firstQuery || (prevQID != curQID))
            {
                firstQuery = false;
//                                // check if qid's are in ascending order
//                                if(prevQID > curQID) 
//                                {
//                                        ostringstream ostr;
//                                        ostr << "QID's of line #" << lineCnt-1 << " and line #" << lineCnt << " are out of order!\n";
//                                        throw CBMRMException(ostr.str(), "CVecFeature::ScanFeatureFile()");
//                                }
            
                numOfAllSubset++;
                subsetSizes->push_back(0);
                prevQID = curQID;
            }
            
            // feature values      
            for(featureCnt=0, featureDimension=0; !iss.eof();)
            {    
                iss >> token;         
                if(IsBlankLine(token)) break;
                if(token[0] == '#') break;
                
                if(sscanf(token.c_str(), scalar_value_format.c_str(), &tmpFval) != 1)
                {
                    ostringstream ostr;
                    ostr << "Feature file contains invalid feature value at line #" << lineCnt+1 << "\n"
                         << "Token: " << token << "\n";
                    throw CBMRMException(ostr.str(), "CVecFeature::ScanFeatureFile()");
                }
                
                // update statistics
                featureDimension++;
                if(SML::abs(tmpFval) > 0.0) 
                    featureCnt++;
            }
            
            // check dimension consistency
            if((prevFeatureDimension != 0) && (prevFeatureDimension != featureDimension)) 
            {
                string msg = "Inconsistent feature vector dimension!\n";
                cout << prevFeatureDimension << "  " << featureDimension << endl;
                throw CBMRMException(msg, "CVecFeature::ScanFeatureFile()");
            }
            
            prevFeatureDimension = featureDimension;
            
            // update statistics
            if(biasFlag)
                featureCnt += 1;  // for shifted hyperplane
            subsetSizes->at(numOfAllSubset-1) += 1;
            nnzSet.push_back(featureCnt);
            allnnz += featureCnt;
            maxNNZ = max(maxNNZ, featureCnt);
            numOfAllExample++;
        }
    }
    
    // for shifted hyperplane: add extra feature (dimension)
    if(biasFlag)
        featureDimension++;
    
    // (whole) data matrix density
    //density = ((double)nnz/featureDimension)/numOfAllExample;
    
    featureFp.close();
}


/**  Load features (the X-part of the dataset) from file into memory.
 */
void CVecFeature::LoadFeatures()
{
    unsigned int* idx = 0;
    double* val = 0;
    unsigned int subsetCnt = 0;
    bool firstQuery = true;
    unsigned int curQID = 0;          // current query id
    unsigned int preQID = 0;         // previous query id
    unsigned int featureCnt = 0;      // number of features per example
    string line = "";
    string token = "";
    ifstream featureFp;
    
    CTimer isstime;
    CTimer addrowtime;
    CTimer creatextime;
    CTimer readlinetime;
    CTimer sscanftime;
    
    // open label file
    featureFp.open(featureFile.c_str());
    
    if(!featureFp.good())
    {    
        string msg = "Cannot open feature vector file <" + featureFile + ">!";
        throw CBMRMException(msg, "CVecFeature::LoadFeatures()");
    }
    
    // allocate memory for data row and subset information
    //idx = new unsigned int[featureDimension];
    //val = new double[featureDimension];
    idx = new unsigned int[maxNNZ];
    val = new double[maxNNZ];
    
    
    // initialize memory content
    memset(idx, 0, sizeof(unsigned int)*maxNNZ);
    memset(val, 0, sizeof(double)*maxNNZ);
    
    creatextime.Start();
    // allocate memory for feature vectors
    if(featureType == DENSE_FEATURE)
        X = new TheMatrix(numOfExample, featureDimension, SML::DENSE);  
    else //if(featureType == SPARSE_FEATURE)
        X = new TheMatrix(numOfExample, featureDimension, SML::SPARSE, nnz);  
    creatextime.Stop();
    
    // skip the examples before #startexample#
    unsigned int skipcnt = 0;
//        do {
//                getline(featureFp,line);
//                trim(line);
//                if(IsBlankLine(line)) continue; // blank line
//                if(line[0] == '#') continue;     // comment                        
//                skipcnt++;
//        } while((skipcnt < startExampleIndex) && !featureFp.eof());
        
        // read feature vectors
    if(featureType == SPARSE_FEATURE)
    {                                  
        // read #nExample# feature vectors
        for(unsigned int lineCnt=0; lineCnt < numOfExample && !featureFp.eof();)
        {      
            // read data line
            readlinetime.Start();
            getline(featureFp, line);
            trim(line);
            readlinetime.Stop();
            if(IsBlankLine(line)) continue;  // blank line
            if(line[0] == '#') continue;     // comment
            
            if(skipcnt < startExampleIndex)
            {
                skipcnt++;
                continue;
            }
            
            isstime.Start();
            istringstream issdata(line);
            isstime.Stop();
            
            // query id is always in front of features !
            size_t pos = line.find("qid:");
            if(pos != string::npos) 
            {
                sscanf(&line[pos], "qid:%u ", &curQID);      
                issdata >> token; // discard this token
            }
            
            if(firstQuery || (preQID != curQID))
            {
                firstQuery = false;
                subsetCnt++;
                preQID = curQID;
                subset[subsetCnt-1].ID = curQID;      
            }                        
            
            //sscanftime.Start();
            // read a feature vector
            for(featureCnt=0; !issdata.eof();)
            {
                issdata >> token;
                if(token[0] == '#') break;
                
                // SPARSE_FEATURE_ELEMENT := "%d:%lf" 
                sscanf(token.c_str(), svec_feature_index_and_value_format.c_str(), &idx[featureCnt], &val[featureCnt]);
                idx[featureCnt] -= 1; // convert feature indexing from 1-based to 0-based
                featureCnt++;
                //printf("(%d,%f)\t",idx[featureCnt-1],val[featureCnt-1]);
            }
            //printf("\n");
            //sscanftime.Stop();
            
            // add extra feature (dimension) for shifted hyperplane
            if(biasFlag)
            {
                idx[featureCnt] = featureDimension-1;  // 0-based indexing
                val[featureCnt] = 1.0;
                featureCnt++;
            }
            
            addrowtime.Start();
            // add row into feature matrix
            X->AddRow(lineCnt, featureCnt, val, idx);
            addrowtime.Stop();
            
            // increase counters
            lineCnt++;
        }
    }
    else  // featureType == DENSE_FEATURE
    { 
        // read #nExample# feature vectors
        for(unsigned int lineCnt=0; lineCnt < numOfExample and !featureFp.eof();)
        {
            // read data line    
            getline(featureFp, line);
            trim(line);
            if(IsBlankLine(line)) continue;  // blank line
            if(line[0] == '#') continue;     // comment
            
            if(skipcnt < startExampleIndex)
            {
                skipcnt++;
                continue;
            }
            
            istringstream issdata(line);
            
            // query id is always in front of features !
            size_t pos = line.find("qid:");
            if(pos != string::npos) 
            {
                sscanf(&line[pos], "qid:%u ", &curQID);      
                issdata >> token; // throw away this token
            }
            
            if(firstQuery || (preQID != curQID))
            {
                firstQuery = false;
                subsetCnt++;
                preQID = curQID;
                subset[subsetCnt-1].ID = curQID;      
            }
            
            
            // read a feature vector
            for(featureCnt=0; !issdata.eof();)
            {
                issdata >> token;
                if(token[0] == '#') break;
                
                // SPARSE_FEATURE_ELEMENT := "%d:%lf" 
                sscanf(token.c_str(), scalar_value_format.c_str(), &val[featureCnt]);
                featureCnt++;
            }
            
            // add extra feature (dimension) for shifted hyperplane
            if(biasFlag) 
            {
                val[featureCnt] = 1.0;
                featureCnt++;
            }
            
            // add row into feature matrix
            X->AddRow(lineCnt, featureCnt, val);
            
            // increace counters
            lineCnt++;
        }
    }
    
    // create row view of X (when requested)   
    if(featureMatrixRowView)
        CreateFeatureMatrixRowViews();
    
    
    // detailed verbosity
    if(vecfeature_verbosity >= 2) 
    {
        double norm = 0;
        X->Norm2(norm);
        cout << "Feature matrix Frobenius norm: " << norm << endl;
    }
    
    if(vecfeature_verbosity >= 2)
    {
        cout << "creatextime : " << creatextime.CPUTotal() << endl;
        cout << "readlinetime : " << readlinetime.CPUTotal() << endl;
        cout << "isstime : " << isstime.CPUTotal() << endl;
        cout << "sscanftime : " << sscanftime.CPUTotal() << endl;
        cout << "addrowtime : " << addrowtime.CPUTotal() << endl;
    }
    
    // clean up
    if(idx) { delete [] idx; idx = 0; }
    if(val) { delete [] val; val = 0; }
    if(subsetSizes) { delete subsetSizes; subsetSizes = 0; }
    
	featureFp.close();
}


/**  Create matrix row view for every feature vector.
 *   These views can be set and get but not shrunk and stretched in actual size
 */
void CVecFeature::CreateFeatureMatrixRowViews()
{
    // feature matrix must be defined first
    // exit if row views had been created before
    if(featureMatrixRowView == true || X == 0) 
        return;
    
    featureMatrixRowView = 1;
    //i.e. x[i] is the i-th feature vector
    x = (TheMatrix**)malloc(sizeof(TheMatrix*)*numOfExample);
    for(unsigned int i=0; i<numOfExample; i++) 
        x[i] = X->CreateMatrixRowView(i);           
}


void CVecFeature::ComputeF(const TheMatrix& w, TheMatrix& f)
{
    X->Dot(w, f); // f = Xw         
}


void CVecFeature::ComputeFi(const TheMatrix& w, TheMatrix& f, const unsigned int i)
{
    assert(featureMatrixRowView);
    assert(i >= 0 && i < numOfExample);
    //x[i]->Dot(w, f); // f_i = x_i*w 
    w.Dot(*(x[i]), f);
}

void CVecFeature::ComputeFi(const TheMatrix& w, double & f, const unsigned int i)
{
    assert(featureMatrixRowView);
    assert(i >= 0 && i < numOfExample);
    w.Dot(*(x[i]), f);
}


void CVecFeature::Grad(const TheMatrix& l, TheMatrix& grad)
{
    //cout<<"X rows: " << (*X).Rows()<<endl;
    //cout<<"X cols: " << (*X).Cols()<<endl;
    //cout<<"l rows: " << l.Rows()<<endl;
    //cout<<"l cols: " << l.Cols()<<endl;
    //cout<<"grad rows: " << grad.Rows()<<endl;
    //cout<<"grad cols: " << grad.Cols()<<endl;
    
    // grad = l*X
    l.Dot(*X, grad); 
    //X->TransposeDot(l, grad); 
    // chteo: X^T*l is slightly faster than l*X in 1-hyperplane cases
    
    //cout<<" printing grad in Grad ... " <<endl;
    //grad.Print();
}


void CVecFeature::AddElement(TheMatrix& w, const unsigned int &i, double scale=1.0)
{
    assert(featureMatrixRowView);
    assert(i >= 0 and i < numOfExample);
    w.ScaleAdd(scale,*(x[i]));
}

#endif
