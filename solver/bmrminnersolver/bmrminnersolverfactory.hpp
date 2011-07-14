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
 *
 * Created: (28/11/2007) 
 *
 * Last Updated:11/07/08
 */

#ifndef _BMRMINNERSOLVERFACTORY_HPP_
#define _BMRMINNERSOLVERFACTORY_HPP_

#include <string>

#include "model.hpp"
#include "bmrminnersolver.hpp"
#include "configuration.hpp"

#include "l2n2_daifletcherpgm.hpp"
#include "l2n2_prloqo.hpp"
#include "l2n2_linesearch.hpp"

#ifdef HAVE_L1N1_INNER_SOLVER
#include "l1n1_clp.hpp"
#endif

/**
 * Factory class for creating new bmrminnersolver object.
 *
 * When you subclass the CBMRMInnerSolver class, you must add an entry into
 * this factory class.
 */
class CBMRMInnerSolverFactory 
{
   public:
      
      /**
       * Return an instance of the bmrminnersolver based on user's argument in the 
       * configuration file.
       *
       */
      static CBMRMInnerSolver* GetBMRMInnerSolver(CModel &model, double lambda)
      {
         CBMRMInnerSolver* innerSolver = 0;
         Configuration &config = Configuration::GetInstance();
         std::string innerSolverType = "";
         
         if(config.IsSet("BMRM.innerSolverType"))
            innerSolverType = config.GetString("BMRM.innerSolverType");
         else
            throw CBMRMException("No BMRM inner solver specified", "CBMRMInnerSolverFactory::GetBMRMInnerSolver()");
         
         // select the innersolver specified by user (in configuration file)
         if(innerSolverType == "L2N2_DaiFletcherPGM")
         {
            innerSolver = new CL2N2_DaiFletcherPGM(lambda);
         }
         else if(innerSolverType == "L2N2_prLOQO")
         {
            innerSolver = new CL2N2_prLOQO(lambda);
         }
         else if(innerSolverType == "L2N2_LineSearch")
         {
            innerSolver = new CL2N2_LineSearch(lambda);
         }
#ifdef HAVE_L1N1_INNER_SOLVER
         else if(innerSolverType == "L1N1_Clp")
         {
            int wLength = model.GetW().Length();
            innerSolver = new CL1N1_Clp(lambda, wLength);
         }
#endif
         else 
         {
            throw CBMRMException("unknown innerSolverType <" + innerSolverType + ">", 
                                 "CBMRMInnerSolverFactory::GetBMRMInnerSolver()");
         }
         
         return innerSolver;
      }
};

#endif
