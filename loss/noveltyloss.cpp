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
 * Created: 18/02/2009
 *
 * Last Updated: 
 */

#include "noveltyloss.hpp"
#include "configuration.hpp"
#include "sml.hpp"

using namespace std;

/** Constructor
 *
 *  @param w [read] weight vector
 *  @param data [read] pointer to dataset
 */
CNoveltyLoss::CNoveltyLoss(CModel* &model, CVecData* &data)
   : CBinaryClassificationLoss(model, data),
     rho(1.0)
{
   // read loss function parameters
   Configuration &config = Configuration::GetInstance();
   
   if(config.IsSet("NOVELTY_DETECTION.rho")) 
      rho = config.GetDouble("NOVELTY_DETECTION.rho");
   
   if(rho <= 0.0)
      throw CBMRMException("Error: NOVELTY_DETECTION.rho must be > 0",
                           "CNoveltyLoss::CNoveltyLoss()");

   if(verbosity > 0)
   {
      cout << "In NOVELTY_DETECTION loss:" << endl;
      cout << "rho = " << rho << endl;
   }
}



/**  
 *  Compute novelty detection loss. 
 *  CAUTION: f is passed by reference and is changed within this function. 
 *  This is done for efficiency reasons, otherwise we would have had to create 
 *  a new copy of f.
 *   
 *  @param loss [write] loss value computed.
 *  @param f [read/write] prediction vector. 
 */
void CNoveltyLoss::Loss(double& loss, TheMatrix& f)
{
   LossAndGrad(loss,f,*l);
}


/**  
 *  Compute loss and gradient of novelty detection loss. 
 *  CAUTION: f is passed by reference and is changed within this
 *  function. This is done for efficiency reasons, otherwise we would
 *  have had to create a new copy of f.
 *   
 *  @param loss [write] loss value computed.
 *  @param f [read/write] prediction vector. 
 *  @param l [write] partial derivative of loss function w.r.t. f
 */
void CNoveltyLoss::LossAndGrad(double& loss, TheMatrix& f, TheMatrix& l)
{
   double* f_array = f.Data();  // pointer to memory location of f (faster element access)
   int len = f.Length();
   l.Zero();  // grad := l'*X
   
   for(int i=0; i < len; i++) 
   {
      if(rho > f_array[i])
      {
         loss += rho - f_array[i];
         l.Set(i, -1.0);
      }
   }
}
