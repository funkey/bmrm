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
 * 			Jin Yu (jin.yu@nicta.com.au)
 *
 * Created: 04/06/2008
 *
 * Last Updated: 19/01/2009
 */

#ifndef _LINESEARCH_MULTILABELLOSS_HPP_
#define _LINESEARCH_MULTILABELLOSS_HPP_

#include "sml.hpp"
#include "multilabelvecdata.hpp"
#include "timer.hpp"
#include "model.hpp"
#include "multilabelloss.hpp"
#include "multilabel_linesearch.hpp"

#define EXACT_LS 0
#define INEXACT_LS 1

/**
 * Class for multilabel classification loss with line-search procedure.
 *
 */

class CLinesearch_MultilabelLoss: public CMultilabelLoss 
{

   protected:
      /** F := [W*x_1; W*x_2; ...; W*x_m] where W is k-by-d and x_i is d-by-1
       */
      TheMatrix *F;
      
      /** i-th matrix row view of F
       */
      TheMatrix *F_i;
      
      /** F := [Wb*x_1; Wb*x_2; ...; Wb*x_m] where Wb is k-by-d and x_i is d-by-1
       */
      TheMatrix *Fb;
      
      /** DF := [P*x_1; P*x_2; ...; P*x_m] where P is k-by-d and x_i is d-by-1
       */
      TheMatrix *DF;
      
      /** w_best, a vector of length k*d
       */
      TheMatrix *Wb;
      
      /** direction P, a vector of length k*d
       */
      TheMatrix *P;
      
      CMultilabel_Linesearch *ls;
      
      /**
       * Regularization constant
       */
      double bmrmLambda;
      
      /**
       * Coefficient used to decide the position of
       * the the new cut:
       * W_{t+1} = theta*W_{t} + (1-theta) Wb
       */
      double theta;
      
      /**
       * Loss evaluated at the best so far iterate
       */
      double loss_of_wb;
      
      /**
       * Step size given by the line search
       */
      double eta;
      
      /**
       * type of the line search method used
       */
      unsigned int ls_method;
      
   public:
      CLinesearch_MultilabelLoss(CModel* &model, CMultilabelVecData* &data);
      
      virtual ~CLinesearch_MultilabelLoss();
      
      // Methods
      /**
       * Calculate loss and gradient
       *
       * @param loss [write] averged sum of losses
       * @param grad [write] a (sub)gradient
       */
      virtual void ComputeLossAndGradient(double& loss, TheMatrix& grad);
      
      /** Compute loss with give weight vector
       *  @param w [read] the given weight vector
       */
      virtual double ComputeLoss(TheMatrix& w);
      
      /**
       * Return loss at the best so far iterate
       */
      virtual double GetLossOfWbest() 
      {
         return loss_of_wb;
      }
      
      /**
       * Return the best so far solution
       */
      virtual void GetWbest(TheMatrix &w_return) 
      {
         w_return.Assign(*Wb);
      }
      
      /**
       * Return the step produced by line search routine
       */
      virtual double GetLineSearchStepSize() 
      {
         return eta;
      }

};

#endif
