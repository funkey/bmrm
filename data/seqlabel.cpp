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
 * Created: (28/01/2008) 
 *
 * Last Updated:
 */

#ifndef _SEQLABEL_CPP_
#define _SEQLABEL_CPP_

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include "common.hpp"
#include "sml.hpp"
#include "configuration.hpp"
#include "bmrmexception.hpp"
#include "timer.hpp"
#include "seqlabel.hpp"

using namespace std;
   
/**  Constructor
 */
CSeqLabel::CSeqLabel()
   : seqlabel_verbosity(0),
     numOfAllLabel(0),
     startExample(0),
     labelFile("")     
{ 
   CTimer loadlabeltime;
   
   // get configurations
   Configuration &config = Configuration::GetInstance();
   
   labelFile = config.GetString("Data.labelFile");
   if(config.IsSet("Data.verbosity"))
      seqlabel_verbosity = config.GetInt("Data.verbosity");
   
   // read labels into memory
   if(seqlabel_verbosity >= 1)
      cout << "Loading label file... "<< endl;
   
   loadlabeltime.Start();
   LoadLabels();
   loadlabeltime.Stop();
   
   // for serial computation, we don't split dataset
   numOfLabel = numOfAllLabel;
   
   if(seqlabel_verbosity >= 2)
      cout << "loadlabeltime   : " << loadlabeltime.CPUTotal() << endl;
}


/**  Load labels (the Y-part of the dataset) from file into memory.
 */
void CSeqLabel::LoadLabels()
{
   string line;
   ifstream labelFp;
   
   // open label file
   labelFp.open(labelFile.c_str());
   
   if(!labelFp.good()) 
   {
      string msg = "Cannot open label file <" + labelFile + ">!";
      throw CBMRMException(msg, "CSeqLabel::LoadLabels()");
   }
   
   numOfAllLabel = 0;
   Y.clear();
   seqlabel_struct tmpseq;    
   
   // read labels  
   while( !labelFp.eof())
   {
      // skip lines with comments
      do {
         getline(labelFp, line);
         trim(line);
      } while((IsBlankLine(line) || line[0] == '#') && !labelFp.eof());
      
      if(labelFp.eof() && IsBlankLine(line)) break;
      
      // clear the temp seq
      tmpseq.ID = 0;
      tmpseq.pos.clear();
      tmpseq.type.clear();
      
      istringstream isslabel(line);
      string token;
      isslabel >> token;   
      
      if(sscanf(token.c_str(),"pid:%u",&(tmpseq.ID)) != 1)
      {
         ostringstream oss;
         oss << "Unable to read pid at line #" << numOfAllLabel+1 << ".";
         throw CBMRMException(oss.str(),"CSeqLabel::LoadLabels()");
      }
      
      
      while(! isslabel.eof())
      {
         isslabel >> token;
         unsigned int tmppos = 0;
         unsigned int tmptype = 0;
         
         if(sscanf(token.c_str(),"%u:%u",&tmppos, &tmptype) == 2)
         {
            tmpseq.pos.push_back(tmppos);
            tmpseq.type.push_back(tmptype);
         }
         else
         {
            ostringstream oss;
            oss << "Unable to read pos:type token at line #" << numOfAllLabel+1 << ".";
            throw CBMRMException(oss.str(),"CSeqLabel::LoadLabels()");
         }                                                                        
      }
      
      // keep only label with at least one sequence element
      if(tmpseq.pos.size() >= 1)
          Y.push_back(tmpseq);
       
      // increace counters
      numOfAllLabel++;
   }
   labelFp.close();
}

void CSeqLabel::Dump()
{
   ofstream ofp("seqlabel.dump");
   for(unsigned int i=0; i<Y.size(); i++)
   {
      ofp << "pid:" << Y[i].ID << " ";
      for(unsigned int j=0; j<Y[i].pos.size(); j++)
      {
         ofp << Y[i].pos[j] << ":" << Y[i].type[j] << " ";
      }
      ofp << endl;
   }
   ofp.close();
}

#endif
