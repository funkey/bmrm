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
 * Created: (14/01/2008) 
 *
 * Last Updated:
 */

#ifndef _ROCSCORELOSS_HPP_
#define _ROCSCORELOSS_HPP_

#include <vector>

#include "common.hpp"
#include "sml.hpp"
#include "binaryclassificationloss.hpp"
#include "model.hpp"

/**  Class to represent binary ROC score loss function:
 */
class CROCScoreLoss : public CBinaryClassificationLoss 
{
   protected:
      /** Use additive label loss?
       *  Used in ramp-loss only
       *  For normal (i.e. convex) ROC Score loss, this flag is always true
       */
      bool _additiveLabelLoss;
      
      /** Number of positive training labels
       */
      unsigned int m_pos;
      
      /** Number of negative training labels
       */
      unsigned int m_neg;                
      
      /** Index vector to be used in sorting of f := X*w
       */
      int* orig_idx;
      int* idx;
      
      void Loss(double& loss, TheMatrix& f);
      void LossAndGrad(double& loss, TheMatrix& f, TheMatrix& l);
      
   public:    
      
      CROCScoreLoss(CModel* &model, CVecData* &data, bool addtiveLabelLoss=true);		
      virtual ~CROCScoreLoss();
};

#endif
