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
 * Last Updated: 15/01/2009
 */

#include "leastsquaresloss.hpp"
#include "sml.hpp"

using namespace std;


/**  
 *  Compute LeastSquare loss. CAUTION: f is passed by reference and is
 *  changed within this function. This is done for efficiency reasons,
 *  otherwise we would have had to create a new copy of f. 
 *   
 *  @param loss [write] loss value computed.
 *  @param f [read/write] prediction vector. 
 */
void CLeastSquaresLoss::Loss(double& loss, TheMatrix& f)
{
   LossAndGrad(loss,f,*l);
}


/**  
 *  Compute loss and partial derivative of LeastSquare loss w.r.t f
 *   
 *  @param loss [write] loss value computed.
 *  @param f [r/w] = X*w
 *  @param l [write] partial derivative of loss w.r.t. f
 */
void CLeastSquaresLoss::LossAndGrad(double& loss, TheMatrix& f, TheMatrix& l)
{
   loss = 0;
   l.Zero();  
   double *Y_array = _data->labels().Data();
   double* f_array = f.Data();
   int len = f.Length();
   for(int i=0; i < len; i++)
   {
      double f_minus_y = (f_array[i] - Y_array[i]);
      loss += f_minus_y * f_minus_y;
      l.Set(i, 2*f_minus_y);
   }
}
