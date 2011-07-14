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
 * Last Updated: (19/11/2007)   
 */

#include "univariateregressionloss.hpp"

using namespace std;

/** Transform function value f := <w,x> into label
 *
 *  @param f [read/write] function values / labels
 */
void CUnivariateRegressionLoss::Transform(TheMatrix& f)
{
   // No transformation
   return;
}


/** Evaluate the performance of the model on _data
 *
 *  @param f [read] function values
 *  @param predict [write] prdicted labels
 */
void CUnivariateRegressionLoss::Perf(TheMatrix& f, TheMatrix& predict)
{
   // Work with a copy of the prediction
   TheMatrix fval(predict);
	
   // compute mse. 
   TheMatrix Y = _data->labels();
   fval.Minus(Y);
   double mse = 0;
   fval.Norm2(mse);
   mse = (mse*mse)/(1.0 * fval.Length());
   
   // dump performanace to stdout
   cout << endl << "Performance:" << endl;
   cout << "1. MSE:  " << mse << endl;
   cout << "2. RMSE: " << sqrt(mse) << endl;
}
