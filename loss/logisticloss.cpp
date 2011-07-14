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
 * Created: (21/10/2008) 
 *
 * Last Updated:
 */

#ifndef _LOGISTICLOSS_CPP_
#define _LOGISTICLOSS_CPP_

#include "logisticloss.hpp"

using namespace std;

const double LN2=0.69314718055994530941;  // natural logarithm


/**  
 *  Compute hinge loss. CAUTION: f is passed by reference and is
 *  changed within this function. This is done for efficiency reasons,
 *  otherwise we would have had to create a new copy of f. 
 *   
 *  @param loss [write] loss value computed.
 *  @param f [read/write] prediction vector. 
 */
void CLogisticLoss::Loss(double& loss, TheMatrix& f)
{
	loss = 0;
	f.ElementWiseMult(_data->labels());  // f = y*f
	double* f_array = f.Data();  // pointer to memory location of f (faster element access)
	int len = f.Length();
	for(int i=0; i < len; i++)
    {
		if(fabs(f_array[i]) == 0.0)
            loss += LN2;
        else if (f_array[i] > 0.0)
            loss += log(1+exp(-f_array[i]));
        else
            loss += log(1+exp(f_array[i])) - f_array[i];
    }
}

/**  
 *  Compute loss and partial derivative of hinge loss w.r.t f
 *   
 *  @param loss [write] loss value computed.
 *  @param f [r/w] = X*w
 *  @param l [write] partial derivative of loss w.r.t. f
 */
void CLogisticLoss::LossAndGrad(double& loss, TheMatrix& f, TheMatrix& l)
{
    l.Zero();  // for gradient computation i.e. grad := l'*X
    f.ElementWiseMult(_data->labels());
    double* f_array = f.Data();  // pointer to memory location of f (faster element access)
    int len = f.Length();	
    double exp_yf = 0.0;

    for(int i=0; i < len; i++)
    {
	if(fabs(f_array[i]) == 0.0)
        {
            loss += LN2;
            l.Set(i,-0.5);
        }
        else if (f_array[i] > 0.0)
        {
            exp_yf = exp(-f_array[i]);
            loss += log(1+exp_yf);
            l.Set(i,-exp_yf/(1+exp_yf));
        }
        else
        {
            exp_yf = exp(f_array[i]);
            loss += log(1+exp_yf) - f_array[i];
            l.Set(i,-1.0/(1+exp_yf));
        }
    }	
    l.ElementWiseMult(_data->labels());
}

#endif
