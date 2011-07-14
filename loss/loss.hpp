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
 * Created: 02/11/2007
 *
 * Last Updated: 16/01/2009
 */

#ifndef _LOSS_HPP_
#define _LOSS_HPP_

#include "sml.hpp"
#include "model.hpp"


/**  
 * Abstract base class for encapsulating a loss object.
 *
 * Sublcass this for arbitrary user-defined losses. When you
 * subclass this class you must also add an entry into CLossFactory. 
 */
class CLoss {
   private:        
      /** Read parameters for CLoss
       */
      void GetParameters();
      
   protected:    
      /** Verbosity level used in loss function
       */
      int verbosity;
      
      /** Scaling factor
       */
      double scalingFactor;
      
      /** Pointer to model object.
       *  This model object contains a weight vector w that
       *  the loss function uses to compute loss and gradient.
       */
      CModel* _model;    
      
      
   public:
      CLoss();
      CLoss(CModel* model, unsigned int numOfW, unsigned int dimOfW, bool biasFeature);
      virtual ~CLoss() {}
      
      
      /** Compute the loss value given the weight vector w in model object  
       * 
       *  @param loss [write] value of the loss.
       */
      virtual void ComputeLoss(double& loss)=0;
      
      
      /** Compute the loss value as well as the gradient of the loss function w.r.t weight 
       *  vector in model object.
       * 
       *  @param loss [write] value of the loss.
       *  @param grad [write] value of the gradient.
       */
      virtual void ComputeLossAndGradient(double& loss, TheMatrix& grad)=0;
      
            
      /** Perform performance evaluation (with test labels) given the model object.
       * 
       *  @param model [read] model object containing weight vector (and possibly other extra parameter)
       */
      virtual void Evaluate(CModel *model)=0;    
      
      
      /** Get the individual loss value computed in previous call to ComputeLossAndGradient()
       *  In general, this method is not useful, so return empty vector by default.
       *  
       *  @param indloss [write] indloss[i] = loss(x_i, y_i, w)
       */
      virtual void GetIndividualLoss(double* &indloss, unsigned int &len) 
      { 
         indloss = 0;
         len = 0;
      }
      
      
      /** Additional methods needed by linesearch loss functions
       */
      virtual double GetLossOfWbest()
      {
         return 0;
      }
      
      virtual void GetWbest(TheMatrix &w_return) {}
      
      virtual double GetScalingFactor() 
      {
         return scalingFactor;
      }
      
};

#endif
