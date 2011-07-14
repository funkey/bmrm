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

#ifndef _REGULARIZER_HPP_
#define _REGULARIZER_HPP_


#include "sml.hpp"
#include "model.hpp"
#include "configuration.hpp"



/** Class for regularizer.
 *  Always vectorize the parameter vector before computing regularization value and 
 *  (sub)gradient.
 */
class CRegularizer
{
   protected:
      int verbosity;
      
   public:
      CRegularizer();
      virtual ~CRegularizer() {}
      
      virtual void ComputeReg(CModel& model, double& regVal) = 0;
      virtual void ComputeRegAndGradient(CModel& model, double& regVal, TheMatrix& regGrad) = 0;
};

#endif
