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
 * Last Updated: (13/11/2007)   
 */

#include "regressionloss.hpp"

using namespace std;

/** Transform function value f := <w,x> into label
 *
 *  @param f [read/write] function values / labels
 */
void CRegressionLoss::Transform(TheMatrix& f)
{
  // do nothing :)
}

/** Evaluate the performance of the model on _data
 *
 *  @param f [read] function values
 *  @param predict [write] prdicted labels
 */
void CRegressionLoss::Perf(TheMatrix& f, TheMatrix& predict)
{
   double* y_array = _data->labels().Data();
   double* f_array = f.Data();
   unsigned int len = f.Length(); 
   double lmse = 0;
   for(unsigned int i = 0; i < len; i++) 
   {
	  double temp = y_array[i] - f_array[i];
	  lmse += temp*temp;
   }
   lmse /= len;

   cout << endl << "Performance:" << endl;
   cout << "1. LMSE:   " << lmse << endl;
}
