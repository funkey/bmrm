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
 * Created: 01/09/2008
 *
 * Last Updated:
 */
      
#ifndef _L1N1_CPP_
#define _L1N1_CPP_

#include "l1n1.hpp"


void CL1N1::ComputeReg(CModel& model, double& reg)
{
   reg = 0;
   TheMatrix &w = model.GetW();
   w.Norm1(reg);
}


/** The subgradient is chosen as sgn(w)
 */
void CL1N1::ComputeRegAndGradient(CModel& model, double& reg, TheMatrix& grad)
{
   reg = 0;
   TheMatrix &w = model.GetW();
   w.Norm1(reg);
   grad.Zero();
   for(int i=0; i<w.Length(); i++)
   {
      double val = 0;
      w.Get(i,val);
      grad.Set(i,SML::sgn(val));
   }
}

#endif
