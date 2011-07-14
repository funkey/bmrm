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

#include "exponentialloss.hpp"
#include "configuration.hpp"
#include "sml.hpp"

using namespace std;


/**  
 *  Compute exponential loss.
 *  CAUTION: f is passed by reference and is changed within this function. 
 *  This is done for efficiency reasons, otherwise we would have had to create 
 *  a new copy of f.
 *   
 *  @param loss [write] loss value computed.
 *  @param f [read/write] prediction vector. 
 */
void CExponentialLoss::Loss(double& loss, TheMatrix& f)
{
   LossAndGrad(loss,f,*l);
}


/**  
 *  Compute loss and gradient of exponential loss. 
 *  CAUTION: f is passed by reference and is changed within this
 *  function. This is done for efficiency reasons, otherwise we would
 *  have had to create a new copy of f.
 *   
 *  @param loss [write] loss value computed.
 *  @param f [read/write] prediction vector. 
 *  @param l [write] partial derivative of loss function w.r.t. f
 */
void CExponentialLoss::LossAndGrad(double& loss, TheMatrix& f, TheMatrix& l)
{
   f.ElementWiseMult(_data->labels());
   f.Scale(-1.0);
   double* yf = f.Data();  // pointer to memory location of f (faster element access)
   double* Y = _data->labels().Data();
   int len = f.Length();
   loss = 0.0;
   l.Zero();  // grad := l'*X
   double max_yf = yf[0];
   for(int i=0; i < len; i++) 
      if(max_yf > yf[i])
         max_yf = yf[i];
   
   for(int i=0; i < len; i++) 
   {
      loss += exp(yf[i]-max_yf);
      l.Set(i, -Y[i]*exp(yf[i]));
   }
   loss *= exp(max_yf);
}
