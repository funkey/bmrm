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
 * Created: (02/11/2007) 
 *
 * Last Updated: (12/01/2009)   
 */

#ifndef _SQUAREDHINGELOSS_CPP_
#define _SQUAREDHINGELOSS_CPP_

#include "squaredhingeloss.hpp"
#include "sml.hpp"

using namespace std;


/**  
 *  Compute squared hinge loss. CAUTION: f is passed by reference and is
 *  changed within this function. This is done for efficiency reasons,
 *  otherwise we would have had to create a new copy of f. 
 *   
 *  @param loss [write] loss value computed
 *  @param f [read/write] prediction vector 
 */
void CSquaredHingeLoss::Loss(double& loss, TheMatrix& f)
{
   loss = 0;
   f.ElementWiseMult(_data->labels());
   double* f_array = f.Data();  // pointer to memory location of f (faster element access)
   int len = f.Length();
   for(int i=0; i < len; i++)
   {
      double val = 1.0 - f_array[i];
      if(val > 0.0) 
         loss += val*val;
   }
   loss *= 0.5;
}


/**  
 *  Compute loss and gradient of squared hinge loss. 
 *   
 *  @param loss [write] loss function computed 
 *  @param f [read] X*w
 *  @param l [write] gradient computed w.r.t. f
 */
void CSquaredHingeLoss::LossAndGrad(double& loss, TheMatrix& f, TheMatrix& l)
{
   loss = 0.0;
   l.Zero();  // for gradient computation i.e. grad := l'*X
   const TheMatrix &Y = _data->labels();
   double* f_array = f.Data();  // pointer to memory location of f (faster element access)
   int len = f.Length();
   double y_i = 0.0;
   
   for(int i=0; i < len; i++) 
   {
      Y.Get(i,y_i);
      double val = 1.0 - y_i*f_array[i];
      if(val > 0.0) 
      {
         loss += val*val;
         l.Set(i, f_array[i] - y_i);
      }
   }
   loss *= 0.5;
}

#endif
