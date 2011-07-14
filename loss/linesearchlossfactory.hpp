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
 * Created: 29/05/2008
 *
 * Last Updated: 16/01/2009
 */

#ifndef _LINESEARCHLOSSFACTORY_HPP_
#define _LINESEARCHLOSSFACTORY_HPP_

#include <string>

#include "common.hpp"
#include "configuration.hpp"
#include "bmrmexception.hpp"
#include "data.hpp"
#include "vecdata.hpp"
#include "multilabelvecdata.hpp"
#include "model.hpp"

// loss function headers
#include "loss.hpp"
#include "linesearch_hingeloss.hpp"
#include "linesearch_multilabelloss.hpp"


/**  Factory class for creating new linesearch Loss instances. 
 */
class CLinesearchLossFactory
{
   public:
      
      /**  Return an instance of loss function based on user's argument in configuration file
       *
       *   @param model [read] Pointer to model object
       *   @param data [read] Pointer to data object
       *   @return loss object
       */
      static CLoss* GetLoss(CModel* &model, CData* &data)
      {         
         CLoss* loss = 0;
         Configuration &config = Configuration::GetInstance();
         
         // select the loss function specified by user (in configuration file)
         if(config.IsSet("Loss.lossFunctionType"))
         {
            std::string lossFunctionType = config.GetString("Loss.lossFunctionType");
            
            if(lossFunctionType == "LINESEARCH_HINGE")
            { 
               CVecData *vecdata = 0;
               if(! (vecdata = dynamic_cast<CVecData*>(data))) 
               {
                  throw CBMRMException("unable to cast data into CVecData",
                                       "CLineSearchLossFactory::GetLoss()");
               }
               loss = new CLinesearch_HingeLoss(model, vecdata);
            }
            else if(lossFunctionType == "LINESEARCH_MULTI_LABEL_CLASSIFICATION")
            { 
               CMultilabelVecData *mlvecdata = 0;
               if(! (mlvecdata = dynamic_cast<CMultilabelVecData*>(data))) 
               {
                  throw CBMRMException("unable to cast data into CVecData",
                                       "CLineSearchLossFactory::GetLoss()");
               }
               loss = new CLinesearch_MultilabelLoss(model, mlvecdata);
            }
            
            else
            {
               throw CBMRMException("ERROR: unrecognised loss function ("+lossFunctionType+")\n", 
                                    "CLineSearchLossFactory::GetLoss()");
            }
         }
         else
         {
            throw CBMRMException("ERROR: no loss function specified!\n", "CLineSearchLossFactory::GetLoss()");
         }
         return loss;  
      }      
};

#endif
