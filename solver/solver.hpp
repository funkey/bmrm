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
 * Created: (25/07/2007) 
 *
 * Last Updated: (06/11/2007)   
 */

#ifndef _SOLVER_HPP_
#define _SOLVER_HPP_

#include "loss.hpp"
#include "model.hpp"
#include "sml.hpp"


/**   
 * Abstract class for solvers which solve the a regularized risk
 * minimization problem. 
 */
class CSolver 
{
   protected:
      CModel* _model;
      CLoss* _loss;
      
   public:      
      CSolver(CModel *model, CLoss *loss) : _model(model), _loss(loss) {}
      virtual ~CSolver() {}
      virtual void Train() = 0;  
};

#endif
