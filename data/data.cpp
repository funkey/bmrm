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
 *          S V N Vishwanathan (SVN.Vishwanathan@nicta.com.au)
 *
 * Created: (02/11/2007) 
 *
 * Last Updated: (06/11/2007)   
 */

#ifndef _DATA_CPP_
#define _DATA_CPP_

#include "data.hpp"
#include "common.hpp"
#include "configuration.hpp"
#include "bmrmexception.hpp"

CData::CData()
   : verbosity(0)
{
        // get configurations
        Configuration &config = Configuration::GetInstance();
   
        if(config.IsSet("Data.verbosity"))
                verbosity = config.GetInt("Data.verbosity");
        else
                verbosity = 0;
}

#endif
