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

#include "poissonloss.hpp"
#include "sml.hpp"

using namespace std;


/**  
 *  Compute Poisson regression loss.
 *  CAUTION: f is passed by reference and is
 *  changed within this function. This is done for efficiency reasons,
 *  otherwise we would have had to create a new copy of f. 
 *   
 *  @param loss [write] loss value computed.
 *  @param f [read/write] prediction vector. 
 */
void CPoissonLoss::Loss(double& loss, TheMatrix& f)
{
   LossAndGrad(loss,f,*l);
}


/**  
 *  Compute loss and gradient of Poisson regression loss w.r.t f
 *   
 *  @param loss [write] loss value computed.
 *  @param f [r/w] = X*w
 *  @param l [write] partial derivative of loss w.r.t. f
 */
void CPoissonLoss::LossAndGrad(double& loss, TheMatrix& f, TheMatrix& l)
{
   loss = 0;
   l.Zero();  
   double *Y_array = _data->labels().Data();
   double* f_array = f.Data();
   int len = f.Length();

   // for(int i=0; i < len; i++)
   // {
   //    loss += exp(f_array[i]) - Y_array[i]*f_array[i];
   //    l.Set(i,exp(f_array[i])-Y_array[i]);
   // }

   double max_f = f_array[0];
   for(int i=1; i < len; i++)
      if(f_array[i] > max_f)
         max_f = f_array[i];

   double part1 = 0.0;
   double part2 = 0.0;
   for(int i=0; i < len; i++)
   {
      part1 += exp(f_array[i]-max_f);
      part2 += Y_array[i]*f_array[i];
      l.Set(i,exp(f_array[i]) - Y_array[i]);
   }
   loss = exp(max_f)*part1 - part2;
}
