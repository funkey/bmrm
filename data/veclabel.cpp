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

#ifndef _VECLABEL_CPP_
#define _VECLABEL_CPP_

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include "common.hpp"
#include "sml.hpp"
#include "configuration.hpp"
#include "bmrmexception.hpp"
#include "timer.hpp"
#include "veclabel.hpp"


using namespace std;
   

CVecLabel::CVecLabel()
   : veclabel_verbosity(0),
     Y(0),
     y(0),
     minLabel(SML::INFTY),
     maxLabel(-SML::INFTY),
     startLabelIndex(0),
     numOfLabel(0),
     numOfAllLabel(0),
     labelType(SINGLE_LABEL),
     labelFile(""),
     labelMatrixRowView(false)
{}



CVecLabel::~CVecLabel()
{  
   if(Y) delete Y;
   if(y) 
   { 
      for(unsigned int i=0; i<numOfLabel; i++)
         delete y[i];
      free(y); 
      y = 0;  
   }
}


/** Load label matrix
 */
void CVecLabel::LoadLabelMatrix(unsigned int startlabel, unsigned int nlabels)
{
   CTimer scanlabeltime;
   CTimer loadlabeltime;
   
   // decide the format string to use
   if(sizeof(double) == sizeof(double))       
      scalar_value_format = "%lf";
   else 
      scalar_value_format = "%f";
   
   // get configurations
   Configuration &config = Configuration::GetInstance();
   
   labelFile = config.GetString("Data.labelFile");
   if(config.IsSet("Data.verbosity"))
      veclabel_verbosity = config.GetInt("Data.verbosity");
   
   // collect some properties of the dataset
   if(veclabel_verbosity >= 1) 
      cout << "Scanning label file... "<< endl;
   
   scanlabeltime.Start();
   ScanLabelFile();
   scanlabeltime.Stop();            
   
   startLabelIndex = startlabel;
   numOfLabel = nlabels;
   
   //printf("startLabelIndex=%d   numOfLabel=%d\n",startLabelIndex,numOfLabel);
   
   // read labels into memory
   if(veclabel_verbosity >= 1)
      cout << "Loading label file... "<< endl;
   
   loadlabeltime.Start();
   LoadLabels();
   loadlabeltime.Stop();
   
   if(veclabel_verbosity >= 2)
   {           
      cout << "scanlabeltime   : " << scanlabeltime.CPUTotal() << endl;
      cout << "loadlabeltime   : " << loadlabeltime.CPUTotal() << endl;
   }        
}


/**  Determine the following properties of the label file:
 *
 *   1.  number of labels
 *   2.  dimensionality of label
 *   3.  largest and smallest label
 *
 *   Note:
 *   1. This method is not quite necessary in serial computation (centralized dataset) mode.
 */
void CVecLabel::ScanLabelFile()
{  
   double tmpLabel = 0;
   unsigned int labelPerEx = 0;
   unsigned int prevLabelPerEx = 0;
   string line = "";
   string token = "";
   ifstream labelFp;
   
   numOfAllLabel = 0;
   
   labelFp.open(labelFile.c_str());
   if(!labelFp.good()) 
   {
      string msg = "Cannot open label file <" + labelFile + ">!";
      throw CBMRMException(msg, "CVecLabel::ScanLabelFile()");
   }
   
   while(!labelFp.eof()) 
   {
      getline(labelFp, line);
      trim(line);
      if(IsBlankLine(line)) continue;  // blank line
      if(line[0] == '#') continue;  // comment line
      istringstream iss(line);
      
      for(labelPerEx = 0; !iss.eof();) 
      {       
         iss >> token;
         if(token[0] == '#') break;
         if(sscanf(token.c_str(), scalar_value_format.c_str(), &tmpLabel) == 1) 
         {
            labelPerEx++;
            minLabel = min(minLabel, tmpLabel);
            maxLabel = max(maxLabel, tmpLabel);
         } 
         else 
         {
            ostringstream ostr;
            ostr << "Label file contains invalid label at line " << numOfAllLabel+1 << "!\n" 
                 << "Token: " << token << "\n";
            throw CBMRMException(ostr.str(), "CVecLabel::ScanLabelFile()");
         }
      }
      
      // ensure all examples have the same number of labels
      assert((prevLabelPerEx == 0) || (prevLabelPerEx == labelPerEx));
      
      prevLabelPerEx = labelPerEx;
      numOfAllLabel++;
   }
   
   if(numOfAllLabel <= 0) 
      throw CBMRMException("Label set is empty!","CVecLabel::ScanLabelFile()");
   
   labelDimension = labelPerEx;
   if(labelDimension > 1) 
      labelType = VECTOR_LABEL;
   
   labelFp.close();
}



/**  Load labels (the Y-part of the dataset) from file into memory.
 */
void CVecLabel::LoadLabels()
{
   unsigned int lineCnt = 0;
   unsigned int labelCnt = 0;        // number of labels per example
   string line;
   ifstream labelFp;
   double tmpLabel;
   
   // open label file
   labelFp.open(labelFile.c_str());
   
   if(!labelFp.good()) 
   {
      string msg = "Cannot open label file <" + labelFile + ">!";
      throw CBMRMException(msg, "CVecLabel::LoadLabels()");
   }
   
   // allocate memory for labelset
   Y = new TheMatrix(numOfLabel, labelDimension, SML::DENSE);
   
   // skip "startLabelIndex" labels
   unsigned int skipcnt = 0;
        
   // read labels  
   for(lineCnt=0; lineCnt < numOfLabel && !labelFp.eof();) 
   {
      // skip lines with comments
      do {
         getline(labelFp, line);
         trim(line);
      } while((IsBlankLine(line) || line[0] == '#') && !labelFp.eof());
      
      if(labelFp.eof() && IsBlankLine(line)) break;
      
      // only read when 
      if(skipcnt < startLabelIndex)
      {
         skipcnt++;
         continue;
      }
      
      istringstream isslabel(line);
      
      if(labelType == SINGLE_LABEL) 
      {
         isslabel >> tmpLabel;         
         Y->Set(lineCnt, tmpLabel);
      }
      else 
      { 
         for(labelCnt=0; !isslabel.eof(); labelCnt++)  
         {
            isslabel >> tmpLabel;
            Y->Set(lineCnt, labelCnt, tmpLabel);
         }
      }
      
      // increace counters
      lineCnt++;
   }
   
   // create row view of Y in the case of VECTOR_LABEL (and when requested)
   if(labelMatrixRowView)
      CreateLabelMatrixRowViews();
   
   labelFp.close();
}


/**  Create matrix row view for every label vector.
 *   These views can be set and get but not shrunk and stretched in actual size
 */
void CVecLabel::CreateLabelMatrixRowViews()
{
   // label must be a vector of length > 1
   // label matrix must be defined first
   // exit if row view had been created before
   if(labelMatrixRowView == true || labelType == SINGLE_LABEL || Y == 0) return;
   
   labelMatrixRowView = 1;
   //i.e. y[i] is the i-th label vector
   y = (TheMatrix**)malloc(sizeof(TheMatrix*)*numOfLabel);
   for(unsigned int i=0; i<numOfLabel; i++)
      y[i] = Y->CreateMatrixRowView(i);  
}

#endif
