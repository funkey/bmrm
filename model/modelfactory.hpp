/* Copyright (c) 2009 NICTA
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
 * Created: (19/12/2007) 
 *
 * Last Updated:
 */

#ifndef _MODELFACTORY_HPP_
#define _MODELFACTORY_HPP_

#include <string>
#include "common.hpp"
#include "configuration.hpp"
#include "bmrmexception.hpp"
#include "model.hpp"


/**  
 * Factory class for creating new Model instances. 
 *
 * When you subclass the CModel class you must also add an entry into
 * this factory class (grep YOUR_MODEL for an example). 
 */
class CModelFactory {
   public:
      
      /**  Return an instance of model based on user's argument in configuration file
       */
      static CModel* GetModel() 
      {         
         CModel* model = 0;
         Configuration &config = Configuration::GetInstance();
         
         // select the loss function specified by user (in configuration file)
         if(config.IsSet("Model.modelType")) 
         {
            std::string modelType = config.GetString("Model.modelType");
            
            if(modelType == "STANDARD")
            {
               model = new CModel();
            }
//          elseif(modelType == "YOUR_MODEL")
//          { 
//             model = new CYourModel();
//          }
//          else
//          {
//             throw CBMRMException("ERROR: unrecognised model ("+modelType+")\n", 
//                                  "CModelFactory::GetModel()");
//          }
            
         } 
         else
         {
            model = new CModel();
         }
         return model;  
      }
      
};

#endif
