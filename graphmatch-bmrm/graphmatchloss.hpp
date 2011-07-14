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
 * Authors: Julian McAucley (julian.mcaculey@nicta.com.au)
 *
 * Created: (27/10/2008) 
 *
 * Last Updated: 17/01/2009
 */

#ifndef _GRAPHMATCHLOSS_HPP_
#define _GRAPHMATCHLOSS_HPP_

#include "common.hpp"
#include "sml.hpp"
#include "data.hpp"
#include "graphdata.hpp"
#include "loss.hpp"
#include <string>
#include "timer.hpp"


/** Class for graphmatch classification loss.
 */
class CGraphMatchLoss : public CLoss
{
   public:
      CGraphMatchLoss(CModel* &model, CGraphData* &data);
      virtual ~CGraphMatchLoss();
      
      // Interfaces
      void ComputeLoss(double& loss);
      void ComputeLossAndGradient(double& loss, TheMatrix& grad);
      void Evaluate(CModel *model);
      double LabelLoss(int* y, int* ybar);
      void LoadModel(std::string modelFilename="");
      void SaveModel(std::string modelFilename="");
      
   protected:
      void Phi(int n, int* y, double* res);
      CGraphData* _data;
      int quad;
      //std::string resultPath;
};

#endif
