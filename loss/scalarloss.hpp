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
 * Created: (02/11/2007) 
 *
 * Last Updated: (13/11/2007)   
 */

#ifndef _SCALARLOSS_HPP_
#define _SCALARLOSS_HPP_

#include <string>

#include "sml.hpp"
#include "model.hpp"
#include "vecdata.hpp"
#include "loss.hpp"

/**  
 * Abstract base class for encapsulating a scalar loss object, that is,
 * a loss which takes as input a real value prediction and a real valued
 * label and outputs a non-negative loss value. Examples include
 * the hinge hinge loss, binary classification loss, and univariate
 *  regression loss.  
 *
 * Sublcass this for arbitrary user-defined scalar losses. When you
 * subclass this class you must also add an entry into CLossFactory and
 * implement the methods Loss and LossAndGrad.
 */

class CScalarLoss : public CLoss 
{
   protected:
      /**  Pointer to a (specific) dataset
       */
      CVecData* _data;
      
      /** Shorthand for number of training examples in the dataset set
       */
      unsigned int m;
      
      /** f := X*w
       */
      TheMatrix* f;
      
      /** l := \partial_{f} lossfunction
       */
      TheMatrix* l;
      
      /** Compute loss function value given f := X*w
       *
       *  @param loss [write] loss function value
       *  @param f [read] := X*w, where X is a matrix of feature vectors and w is weight vector
       */
      virtual void Loss(double& loss, TheMatrix& f)=0;
      
      
      /** Compute loss function value and partial derivative of loss function w.r.t. f, given f := X*w
       *
       *  @param loss [write] loss function value
       *  @param f [read] := X*w, where X is a matrix of feature vectors and w is weight vector
       *  @param l [write] partial derivative of loss function w.r.t. f
       */
      virtual void LossAndGrad(double& loss, TheMatrix& f, TheMatrix& l)=0;
      
      
      /** Transform f into corresponding labels
       *  
       *  @param f [r/w] real values --> labels
       */
      virtual void Transform(TheMatrix& f) = 0;
      
      
      /** Performance measures
       *
       *  @param f [read] X*w
       *  @param predict [write] predicted labels
       */
      virtual void Perf(TheMatrix& f, TheMatrix& predict) = 0;
      
	
   public:      
      CScalarLoss(CModel* &model, CVecData* &data);
      virtual ~CScalarLoss();
      
      virtual void ComputeLoss(double& loss);
      virtual void ComputeLossAndGradient(double& loss, TheMatrix& grad);
      virtual void Evaluate(CModel *model);
};

#endif
