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
 * Created: (26/01/2008) 
 *
 * Last Updated:
 */

#ifndef _MULTILABELLOSS_HPP_
#define _MULTILABELLOSS_HPP_


#include "sml.hpp"
#include "data.hpp"
#include "multilabelvecdata.hpp"
#include "loss.hpp"
#include "model.hpp"

/** Class for max-margin multilabel loss.
 *  
 *  The loss is formally defined as 
 *  loss = max(0, f(x, y*) - f(x, y)) 
 *  where f(x, y') = <w_{y'}, x>, 
 *        y* := argmax_{y'} f(x, y'), y' \in Y' (i.e. all valid y except Y)
 *        y is the true label set Y. 
 */
class CMultilabelLoss : public CLoss 
{    
   protected:
      /** Pointer to data
       */
      CMultilabelVecData* _data;
      
      /** Number of examples in data
       */
      unsigned int m;
      
      /** Number of classes in data
       */
      unsigned int numOfClass;
      
      /** Prediction f := Xw
       */
      TheMatrix *f;
      
      /** Weight vector in matrix form
       */
      TheMatrix *matW;
      
      /** Gradient of loss w.r.t. w in matrix form
       */
      TheMatrix *matG;  // local copy of weight vector and gradient vector in matrix form
      
      /** Explicit row view of matG
       */
      TheMatrix **g;
      
      /** The set of false labels
       */
      std::vector<bool> Ybar;
      
      /** Compute decision function values and predict labels
       *
       *  @param model [read] Model object containing a weight vector
       *  @param ybar [write] Predicted labels
       *  @param f_ybar [write] Decision function values
       */
      virtual void DecisionAndPrediction(CModel* model, int* ybar, double* f_ybar);
      
   public:
      CMultilabelLoss() {}
      CMultilabelLoss(CModel* &model, CMultilabelVecData* &data);
      virtual ~CMultilabelLoss();
      
      // Methods
      virtual void ComputeLoss(double& loss);
      virtual void ComputeLossAndGradient(double& loss, TheMatrix& grad);
      virtual void Evaluate(CModel* model);
};

#endif
