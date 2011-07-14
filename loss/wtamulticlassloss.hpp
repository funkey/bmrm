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
 * Last Updated: (12/01/2009)   
 */

#ifndef _WTAMULTICLASSLOSS_HPP_
#define _WTAMULTICLASSLOSS_HPP_

#include <string>
#include "sml.hpp"
#include "data.hpp"
#include "vecdata.hpp"
#include "loss.hpp"
#include "model.hpp"

/**  
 * Class for winner-takes-all multiclass losses. 
 * 
 * The assumption is that the label set is a small finite set (say 10 in
 * the case of digit recognition), and the loss measures the discrepancy
 * between the largest function value and function value at the true
 * label (see eg. Crammer and Singer 2001 for details). 
 * 
 * The loss is formally defined as 
 * loss = max(0, margin + f(x, y*) - f(x, y) 
 * where f(x, y') = <w_{y'}, x>, y* := argmax_{y' \neq y} f(x, y') and y is the
 * true label. 
 */
class CWTAMulticlassLoss : public CLoss 
{
   protected:
      /** The type of the label loss function 
       */
      std::string labelLossType;
      
      /** The type of the margin scaling function
       */
      std::string marginScalingFunctionType;
    
      /** Pointer to data
       */
      CVecData* _data;
      
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
      
      /** Whether to add label loss to the loss function or not
       *  if true:  loss = max_y' <w, \phi(x,y') - \phi(x,y)> + \Delta(y,y')
       *  if false: loss = max_y' <w, \phi(x,y') - \phi(x,y)>
       */
      bool _additiveLabelLoss;
      
      
      /** array to keep individual loss
       */
      double* _indloss;
      
      
      /** Function pointer to margin scaling function.
       *
       *  @param delta [read] the cost for wrong label prediction
       *  @return margin scaling factor
       */
      double (CWTAMulticlassLoss::*Gamma)(const double& labelloss) const;
      
      
      /** Margin scaling function: identity / no scaling / "slack re-scaling"
       */
      double MarginRescaling(const double& labelloss) const
      {
         return 1;
      }
      
      
      /** Margin scaling function: "margin re-scaling"    
       *  (with 0/1 loss, slack re-scaling is the same as margin re-scaling)
       */
      double SlackRescaling(const double& labelloss) const
      {
         return labelloss;
      }
      
      
      /** Function pointer to label loss function
       *
       *  @param y [read] actual label
       *  @param ybar [read] predicted label
       *  @return label loss 
       */
      double (CWTAMulticlassLoss::*Delta)(const int& y, const int& ybar) const;
      
      
      /** Cost for predicting wrong label: 0/1 loss 
       */
      double ZeroOne(const int& y, const int& ybar) const
      {
         return ((y != ybar) ? 1 : 0);
      }
      
      
      /** Cost for predicting wrong label: 0/100 loss 
       */
      double ZeroOneHundred(const int& y, const int& ybar) const
      {
         return ((y != ybar) ? 100 : 0);
      }
      
      
      /** Cost for predicting wrong label: squared difference 
       */
      double SquaredDifference(const int& y, const int& ybar) const
      {
         return ((y-ybar)*(y-ybar));
      }   
      
      
      /** Compute decision function values and predict labels
       *
       *  @param model [read] Model object containing a weight vector
       *  @param ybar [write] Predicted labels
       *  @param f_ybar [write] Decision function values
       */
      virtual void DecisionAndPrediction(CModel* model, int* ybar, double* f_ybar);
      
   public:
      CWTAMulticlassLoss() {}
      CWTAMulticlassLoss(CModel* &model, CVecData* &data, bool useAdditiveLabelLoss=true);
      virtual ~CWTAMulticlassLoss();
      
      // Methods
      virtual void ComputeLoss(double& loss);
      virtual void ComputeLossAndGradient(double& loss, TheMatrix& grad);
      virtual void Evaluate(CModel* model);
      virtual void GetIndividualLoss(double* &indloss, unsigned int &len);
};

#endif
