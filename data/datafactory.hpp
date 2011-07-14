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
 * Authors: S V N Vishwanathan (SVN.Vishwanathan@nicta.com.au)
 *          Choon Hui Teo (ChoonHui.Teo@anu.edu.au)
 *
 * Created: (02/11/2007) 
 *
 * Last Updated: (15/04/2008)   
 */

#ifndef _DATAFACTORY_HPP_
#define _DATAFACTORY_HPP_

#include <string>

#include "bmrmexception.hpp"
#include "configuration.hpp"
#include "data.hpp"
#include "vecdata.hpp"
#include "multilabelvecdata.hpp"

/**  
 * Factory class for creating new Data instances. 
 *
 * When you subclass the CData class you must also add an entry into
 * this factory class (grep YOUR_DATA_FORMAT for an example). 
 */
class CDataFactory
{
   public:
      /**  Return a data object based on user's argument in configuration file.
       *   For serial computation, start and nparts are dummies.
       *   For distributed/parallel computation,
       *     load a portion of whole dataset: divide dataset into "nparts" parts and 
       *     load only the "start"-th part.
       *
       *   @param start [read] The part of dataset (divided into nparts) this machine should load
       *   @param nparts [read] The number of parts the original dataset will be divided into
       *   @return data object
       */
      static CData* GetData(unsigned int start=0, unsigned int nparts=1)
      {  
         CData *ds = 0;
         Configuration &config = Configuration::GetInstance();
         
         // default to this format
         std::string dataFormat = "VECTOR_LABEL_VECTOR_FEATURE";
         
         // unless the user specifies otherwise in the config 
         if(config.IsSet("Data.format"))
            dataFormat = config.GetString("Data.format");
         
         if(dataFormat == "VECTOR_LABEL_VECTOR_FEATURE")
         {
            ds = new CVecData(start,nparts);            
         }
         else if(dataFormat == "VARIABLE_LENGTH_VECTOR_LABEL_VECTOR_FEATURE")
         {
            ds = new CMultilabelVecData(start,nparts);            
         }
         // else if(dataFormat == "YOUR_DATA_FORMAT")
         //{
         //   ds = new CYourDataFormat();
         //}
         else
         {
            throw CBMRMException("ERROR: unrecognised data format ("+dataFormat+")\n", "CDataFactory::GetData()");
         }
         
         return ds;
      }   
};

#endif
