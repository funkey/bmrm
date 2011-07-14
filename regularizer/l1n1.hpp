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

#ifndef _L1N1_HPP_
#define _L1N1_HPP_


#include "regularizer.hpp"


/** Class for L1 norm regularizer.
 *  Always vectorize the parameter vector before computing regularization value and 
 *  (sub)gradient.
 */
class CL1N1 : public CRegularizer
{
   public:
      CL1N1() : CRegularizer() {}
      virtual ~CL1N1() {}
      
      virtual void ComputeReg(CModel& model, double& regVal);
      virtual void ComputeRegAndGradient(CModel& model, double& regVal, TheMatrix& regGrad);
};

#endif
