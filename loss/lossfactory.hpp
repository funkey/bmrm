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
 * Last Updated: 15/01/2009
 */

#ifndef _LOSSFACTORY_HPP_
#define _LOSSFACTORY_HPP_

#include <string>

#include "common.hpp"
#include "configuration.hpp"
#include "bmrmexception.hpp"
#include "data.hpp"
#include "vecdata.hpp"
#include "model.hpp"

// binary classification
#include "loss.hpp"
#include "hingeloss.hpp"
#include "squaredhingeloss.hpp"
#include "huberhingeloss.hpp"
#include "logisticloss.hpp"
#include "exponentialloss.hpp"

#ifndef PARALLEL_BMRM
#include "rocscoreloss.hpp"
#include "fbetaloss.hpp"
#endif


// univariate regression
#include "epsiloninsensitiveloss.hpp"
#include "leastsquaresloss.hpp"
#include "leastabsdevloss.hpp"
#include "quantileloss.hpp"
#include "poissonloss.hpp"
#include "huberrobustloss.hpp"

// novelty detection
#include "noveltyloss.hpp"


// multiclass/multilabel classification
#include "wtamulticlassloss.hpp"
#include "multilabelloss.hpp"

// ranking
#include "ndcgrankloss.hpp"


/**  
 * Factory class for creating new Loss instances. 
 *
 * When you subclass the CLoss class you must also add an entry into
 * this factory class (grep YOUR_LOSS_FUNCTION for an example). 
 * (dynamic) cast the generic data object handler to your specific data object type, if
 * you have a special subclass of CData.
 */
class CLossFactory
{
   public:
      
      /**  Return an instance of loss function based on user's argument in configuration file
       *
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
            
            if(lossFunctionType == "HINGE")
            { 
               CVecData *vecdata = 0;
               if(! (vecdata = dynamic_cast<CVecData*>(data)))
               {
                  throw CBMRMException("unable to cast data into CVecData",
                                       "CLossFactory::GetLoss()");
               }
               loss = new CHingeLoss(model, vecdata);
            } 
            else if(lossFunctionType == "SQUARED_HINGE")
            { 
               CVecData *vecdata = 0;
               if(! (vecdata = dynamic_cast<CVecData*>(data))) 
               {
                  throw CBMRMException("unable to cast data into CVecData",
                                       "CLossFactory::GetLoss()");
               }
               loss = new CSquaredHingeLoss(model, vecdata);
            }
            else if(lossFunctionType == "HUBER_HINGE")
            { 
               CVecData *vecdata = 0;
               if(! (vecdata = dynamic_cast<CVecData*>(data))) 
               {
                  throw CBMRMException("unable to cast data into CVecData",
                                       "CLossFactory::GetLoss()");
               }
               loss = new CHuberHingeLoss(model, vecdata);
            }
            else if(lossFunctionType == "LOGISTIC")
            { 
               CVecData *vecdata = 0;
               if(! (vecdata = dynamic_cast<CVecData*>(data))) 
               {
                  throw CBMRMException("unable to cast data into CVecData",
                                       "CLossFactory::GetLoss()");
               }
               loss = new CLogisticLoss(model, vecdata);
            }
            else if(lossFunctionType == "EXPONENTIAL")
            { 
               CVecData *vecdata = 0;
               if(! (vecdata = dynamic_cast<CVecData*>(data))) 
               {
                  throw CBMRMException("unable to cast data into CVecData",
                                       "CLossFactory::GetLoss()");
               }
               loss = new CExponentialLoss(model, vecdata);
            }
#ifndef PARALLEL_BMRM
            else if(lossFunctionType == "ROC_SCORE")
            { 
               CVecData *vecdata = 0;
               if(! (vecdata = dynamic_cast<CVecData*>(data))) 
               {
                  throw CBMRMException("unable to cast data into CVecData",
                                       "CLossFactory::GetLoss()");
               }
               loss = new CROCScoreLoss(model, vecdata);
            }
            else if(lossFunctionType == "F_BETA")
            { 
               CVecData *vecdata = 0;
               if(! (vecdata = dynamic_cast<CVecData*>(data))) 
               {
                  throw CBMRMException("unable to cast data into CVecData",
                                       "CLossFactory::GetLoss()");
               }
               loss = new CFBetaLoss(model, vecdata);
            }
#endif
            else if(lossFunctionType == "EPSILON_INSENSITIVE")
            { 
               CVecData *vecdata = 0;
               if(! (vecdata = dynamic_cast<CVecData*>(data)))
               {
                  throw CBMRMException("unable to cast data into CVecData",
                                       "CLossFactory::GetLoss()");
               }
               loss = new CEpsilonInsensitiveLoss(model, vecdata);
            }
            else if(lossFunctionType == "LEAST_SQUARES")
            {
               CVecData *vecdata = 0;
               if(! (vecdata = dynamic_cast<CVecData*>(data))) 
               {
                  throw CBMRMException("unable to cast data into CVecData",
                                       "CLossFactory::GetLoss()");
               }
               loss = new CLeastSquaresLoss(model, vecdata);
            } 
            else if(lossFunctionType == "LEAST_ABSOLUTE_DEVIATION")
            {
               CVecData *vecdata = 0;
               if(! (vecdata = dynamic_cast<CVecData*>(data))) 
               {
                  throw CBMRMException("unable to cast data into CVecData",
                                       "CLossFactory::GetLoss()");
               }
               loss = new CLeastAbsDevLoss(model, vecdata);
            } 
            else if(lossFunctionType == "QUANTILE_REGRESSION")
            {
               CVecData *vecdata = 0;
               if(! (vecdata = dynamic_cast<CVecData*>(data))) 
               {
                  throw CBMRMException("unable to cast data into CVecData",
                                       "CLossFactory::GetLoss()");
               }
               loss = new CQuantileLoss(model, vecdata);
            } 
            else if(lossFunctionType == "POISSON_REGRESSION")
            {
               CVecData *vecdata = 0;
               if(! (vecdata = dynamic_cast<CVecData*>(data))) 
               {
                  throw CBMRMException("unable to cast data into CVecData",
                                       "CLossFactory::GetLoss()");
               }
               loss = new CPoissonLoss(model, vecdata);
            } 
            else if(lossFunctionType == "HUBER_ROBUST_REGRESSION")
            {
               CVecData *vecdata = 0;
               if(! (vecdata = dynamic_cast<CVecData*>(data))) 
               {
                  throw CBMRMException("unable to cast data into CVecData",
                                       "CLossFactory::GetLoss()");
               }
               loss = new CHuberRobustLoss(model, vecdata);
            } 
            else if(lossFunctionType == "NOVELTY_DETECTION")
            {
               CVecData *vecdata = 0;
               if(! (vecdata = dynamic_cast<CVecData*>(data))) 
               {
                  throw CBMRMException("unable to cast data into CVecData",
                                       "CLossFactory::GetLoss()");
               }
               loss = new CNoveltyLoss(model, vecdata);
            } 
            else if(lossFunctionType == "WTA_MULTICLASS")
            {
               CVecData *vecdata = 0;
               if(! (vecdata = dynamic_cast<CVecData*>(data))) 
               {
                  throw CBMRMException("unable to cast data into CVecData",
                                       "CLossFactory::GetLoss()");
               }
               loss = new CWTAMulticlassLoss(model, vecdata);
            }
            else if(lossFunctionType == "MULTI_LABEL_CLASSIFICATION")
            {
               CMultilabelVecData *mlvecdata = 0;
               if(! (mlvecdata = dynamic_cast<CMultilabelVecData*>(data))) 
               {
                  throw CBMRMException("unable to cast data into CMultilabelVecData",
                                       "CLossFactory::GetLoss()");
               }
               loss = new CMultilabelLoss(model, mlvecdata);
            }
            else if(lossFunctionType == "NDCG_RANK")
            {
               CVecData *vecdata = 0;
               if(! (vecdata = dynamic_cast<CVecData*>(data))) 
               {
                  throw CBMRMException("unable to cast data into CVecData",
                                       "CLossFactory::GetLoss()");
               }
               loss = new CNDCGRankLoss(model, vecdata);
            } 
            else
            {
               throw CBMRMException("ERROR: unrecognised loss function ("+lossFunctionType+")\n", 
                                    "CLossFactory::GetLoss()");
            }
         }
         else
         {
            throw CBMRMException("ERROR: no loss function specified!\n", "CLossFactory::GetLoss()");
         }
         return loss;  
      }      
};

#endif
