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

#include "quantileloss.hpp"
#include "configuration.hpp"
#include "sml.hpp"

using namespace std;

/** Constructor
 *
 *  @param w [read] weight vector
 *  @param data [read] pointer to dataset
 */
CQuantileLoss::CQuantileLoss(CModel* &model, CVecData* &data)
   : CUnivariateRegressionLoss(model, data),
     tau(0.5)
{
   // read loss function parameters
   Configuration &config = Configuration::GetInstance();
   
   if(config.IsSet("QUANTILE_REGRESSION.tau")) 
      tau = config.GetDouble("QUANTILE_REGRESSION.tau");
   
   if(tau <= 0.0 || tau >= 1.0)
      throw CBMRMException("Error: QUANTILE_REGRESSION.tau must be in the range (0,1)",
                           "CQuantileLoss::CQuantileLoss()");

   if(verbosity > 0)
   {
      cout << "In QUANTILE_REGRESSION loss:" << endl;
      cout << "tau = " << tau << endl;
   }
}



/**  
 *  Compute quantile regression loss. 
 *  CAUTION: f is passed by reference and is changed within this function. 
 *  This is done for efficiency reasons, otherwise we would have had to create 
 *  a new copy of f.
 *   
 *  @param loss [write] loss value computed.
 *  @param f [read/write] prediction vector. 
 */
void CQuantileLoss::Loss(double& loss, TheMatrix& f)
{
   LossAndGrad(loss,f,*l);
}


/**  
 *  Compute loss and gradient of quantile regression loss. 
 *  CAUTION: f is passed by reference and is changed within this
 *  function. This is done for efficiency reasons, otherwise we would
 *  have had to create a new copy of f.
 *   
 *  @param loss [write] loss value computed.
 *  @param f [read/write] prediction vector. 
 *  @param l [write] partial derivative of loss function w.r.t. f
 */
void CQuantileLoss::LossAndGrad(double& loss, TheMatrix& f, TheMatrix& l)
{
   double* f_array = f.Data();  // pointer to memory location of f (faster element access)
   double* y_array = _data->labels().Data();
   int len = f.Length();
   
   l.Zero();  // grad := l'*X
   
   for(int i=0; i < len; i++) 

   {
      if(f_array[i] >= y_array[i])
      {
         loss += tau*(f_array[i]-y_array[i]);
         l.Set(i, tau);
      }
      else
      {
         loss += (1-tau)*(y_array[i]-f_array[i]);
         l.Set(i, tau-1);
      }
   }
}
