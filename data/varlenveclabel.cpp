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

#ifndef _VARLENVECLABEL_CPP_
#define _VARLENVECLABEL_CPP_

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include "common.hpp"
#include "sml.hpp"
#include "configuration.hpp"
#include "bmrmexception.hpp"
#include "timer.hpp"
#include "varlenveclabel.hpp"


using namespace std;

CVarLenVecLabel::CVarLenVecLabel()
   : varlenveclabel_verbosity(0),
     Y(0),
     avgLabelDimension(0.0),
     minLabel(SML::INFTY),
     maxLabel(-SML::INFTY),
     numOfLabel(0),
     numOfAllLabel(0),
     startLabelIndex(0),
     labelFile(""),
     labelIncrement(0.0)
     
{}


CVarLenVecLabel::~CVarLenVecLabel()
{}


/** Load labels
 */
void CVarLenVecLabel::LoadVarLenLabel(unsigned int startlabel, unsigned int nlabels)
{
   CTimer scanlabeltime;
   CTimer loadlabeltime;
   
   scalar_value_format = "%lf";
   
   // get configurations
   Configuration &config = Configuration::GetInstance();
   
   labelFile = config.GetString("Data.labelFile");
   if(config.IsSet("Data.verbosity"))
      varlenveclabel_verbosity = config.GetInt("Data.verbosity");

   if(config.IsSet("Data.labelIncreasedBy"))
      labelIncrement = config.GetDouble("Data.labelIncreasedBy");
   
   // collect some properties of the dataset
   if(varlenveclabel_verbosity >= 1) 
      cout << "Scanning label file... "<< endl;
   
   scanlabeltime.Start();
   ScanLabelFile();
   scanlabeltime.Stop();            
   
   numOfLabel = nlabels;
   startLabelIndex = startlabel;
   
   // read labels into memory
   if(varlenveclabel_verbosity >= 1)
      cout << "Loading label file... "<< endl;
   
   loadlabeltime.Start();
   LoadLabels();
   loadlabeltime.Stop();
   
   if(varlenveclabel_verbosity >= 2)
   {           
      cout << "scanlabeltime   : " << scanlabeltime.CPUTotal() << endl;
      cout << "loadlabeltime   : " << loadlabeltime.CPUTotal() << endl;
   }        
}


/**  Determine the following properties of the label file:
 *   1.  type of label
 *   2.  number of labels
 *   3.  dimensionality of label
 *   4.  largest label
 */
void CVarLenVecLabel::ScanLabelFile()
{  
   double tmpLabel = 0;
   string line = "";
   string token = "";
   ifstream labelFp;
   
   numOfAllLabel = 0;
   
   labelFp.open(labelFile.c_str());
   if(!labelFp.good()) 
   {
      string msg = "Cannot open label file <" + labelFile + ">!";
      throw CBMRMException(msg, "CVarLenVecLabel::ScanLabelFile()");
   }
   
   while(!labelFp.eof()) 
   {
      getline(labelFp, line);
      trim(line);
      if(IsBlankLine(line)) continue;  // blank line
      if(line[0] == '#') continue;  // comment line
      istringstream iss(line);
      
      //while(!iss.eof()) 
      while(getline(iss,token,',')) 
      {       
         //iss >> token;
         if(token[0] == '#') break;
         if(sscanf(token.c_str(), scalar_value_format.c_str(), &tmpLabel) == 1) 
         {
            tmpLabel += labelIncrement;
            minLabel = min(minLabel, tmpLabel);
            maxLabel = max(maxLabel, tmpLabel);
         } 
         else 
         {
            ostringstream ostr;
            ostr << "Label file contains invalid label at line " << numOfAllLabel+1 << "!\n" 
                 << "Token: " << token << "\n";
            throw CBMRMException(ostr.str(), "CVarLenVecLabel::ScanLabelFile()");
         }
      }     
      numOfAllLabel++;
   }
   
   if(numOfAllLabel <= 0) 
      throw CBMRMException("Label set is empty!","CVarLenVecLabel::ScanLabelFile()");
   
   labelFp.close();
}



/**  Load labels (the Y-part of the dataset) from file into memory.
 */
void CVarLenVecLabel::LoadLabels()
{
   unsigned int lineCnt = 0;
   unsigned int labelCnt = 0;        // number of labels per example
   string line, token;
   ifstream labelFp;
   double tmpLabel;
   vector<double> tmpLabelset;
   
   // open label file
   labelFp.open(labelFile.c_str());
   
   if(!labelFp.good()) 
   {
      string msg = "Cannot open label file <" + labelFile + ">!";
      throw CBMRMException(msg, "CVarLenVecLabel::LoadLabels()");
   }
   
   // read labels  
   unsigned int skipcnt = 0;
   for(lineCnt=0; lineCnt < numOfLabel && !labelFp.eof() ;) 
   {
      // skip lines with comments
      do {
         getline(labelFp, line);
         trim(line);
      } while((IsBlankLine(line) || line[0] == '#') && !labelFp.eof());
      
      if(labelFp.eof() && IsBlankLine(line)) break;
      
      if(skipcnt < startLabelIndex)
      {
         skipcnt++;
         continue;
      }
      
      istringstream isslabel(line);
      
      //for(labelCnt=0; !isslabel.eof(); labelCnt++)  
      labelCnt = 0;
      while(getline(isslabel,token,','))
      {
         //isslabel >> tmpLabel;
         stringstream(token) >> tmpLabel;
         tmpLabel += labelIncrement;
         tmpLabelset.push_back(tmpLabel);
         avgLabelDimension += 1;
         labelCnt++;
      }
      assert(tmpLabelset.size() != 0);
      Y.push_back(tmpLabelset);
      tmpLabelset.clear();
      lineCnt++;
   }
   assert((unsigned int)Y.size() == numOfLabel);
   
   avgLabelDimension /= numOfAllLabel;
   labelFp.close();
}


/** show labels
 */
void CVarLenVecLabel::ShowLabel()
{
   for(unsigned int i=0; i<Y.size(); i++)
   {
      for(unsigned int j=0; j<Y[i].size(); j++)
         printf("%f\t",Y[i][j]);
      printf("\n");
   }
}

#endif
