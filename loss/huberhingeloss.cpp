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

#include "huberhingeloss.hpp"
#include "configuration.hpp"
#include "sml.hpp"

using namespace std;

/** Constructor
 *
 *  @param w [read] weight vector
 *  @param data [read] pointer to dataset
 */
CHuberHingeLoss::CHuberHingeLoss(CModel* &model, CVecData* &data)
   : CBinaryClassificationLoss(model, data),
     h(0.01)
{
   // read loss function parameters
   Configuration &config = Configuration::GetInstance();
   
   if(config.IsSet("HUBER_HINGE.h")) 
      h = config.GetDouble("HUBER_HINGE.h");
   
   if(h <= 0.0)
      throw CBMRMException("Error: HUBER_HINGE.h must be > 0",
                           "CHuberHingeLoss::CHuberHingeLoss()");

   if(verbosity > 0)
   {
      cout << "In HUBER_HINGE loss:" << endl;
      cout << "h = " << h << endl;
   }
}



/**  
 *  Compute Huber hinge loss. 
 *  CAUTION: f is passed by reference and is changed within this function. 
 *  This is done for efficiency reasons, otherwise we would have had to create 
 *  a new copy of f.
 *   
 *  @param loss [write] loss value computed.
 *  @param f [read/write] prediction vector. 
 */
void CHuberHingeLoss::Loss(double& loss, TheMatrix& f)
{
   LossAndGrad(loss,f,*l);
}


/**  
 *  Compute loss and gradient of Huber hinge loss. 
 *  CAUTION: f is passed by reference and is changed within this
 *  function. This is done for efficiency reasons, otherwise we would
 *  have had to create a new copy of f.
 *   
 *  @param loss [write] loss value computed.
 *  @param f [read/write] prediction vector. 
 *  @param l [write] partial derivative of loss function w.r.t. f
 */
void CHuberHingeLoss::LossAndGrad(double& loss, TheMatrix& f, TheMatrix& l)
{
   f.ElementWiseMult(_data->labels());
   double* yf = f.Data();
   double* Y = _data->labels().Data();
   int len = f.Length();
   loss = 0.0;
   l.Zero();

   for(int i=0; i < len; i++) 
   {
      double v = 1-yf[i];
      if(h < v)
      {
         loss += v;
         l.Set(i,-Y[i]);
      }
      else if(-h > v) {}
      else
      {
         loss += (v+h)*(v+h)/4/h;
         l.Set(i, -Y[i]*(v+h)/2/h);
      }
   }
}
