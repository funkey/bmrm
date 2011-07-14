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

#ifndef _REGULARIZERFACTORY_HPP_
#define _REGULARIZERFACTORY_HPP_

#include <string>
#include "configuration.hpp"
#include "regularizer.hpp"
#include "l2n2.hpp"
#include "l1n1.hpp"

/** Factory for creating new regularizer object.
 */
class CRegularizerFactory
{
   public:
      
      static CRegularizer* GetRegularizer()
      {
         CRegularizer *reg = 0;
         Configuration &config = Configuration::GetInstance();
         
         if(config.IsSet("Regularizer.regularizerType"))
         {
			std::string regularizerType = config.GetString("Regularizer.regularizerType");
			if(regularizerType == "L2N2")
			{
               reg = new CL2N2();
			}
			else if(regularizerType == "L1N1")
			{
               reg = new CL1N1();
			}
         }
         else
         {
			reg = new CL2N2();
         }
         return reg;
      }
};

#endif
