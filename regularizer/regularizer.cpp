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

#ifndef _REGULARIZER_CPP_
#define _REGULARIZER_CPP_


#include "regularizer.hpp"
#include "configuration.hpp"

CRegularizer::CRegularizer()
   : verbosity(0)
{
	Configuration &config = Configuration::GetInstance();
    
	if(config.IsSet("Regularizer.verbosity"))
       verbosity = config.GetInt("Regularizer.verbosity");
}

#endif
