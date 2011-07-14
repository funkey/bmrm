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

#ifndef _SCALARLOSS_CPP_
#define _SCALARLOSS_CPP_

#include <fstream>
#include <sstream>

#include "scalarloss.hpp"
#include "configuration.hpp"

using namespace std;

/**  
 * Constructor
 *
 *   @param model [read/write] pointer to model object
 *   @param data [read/write] pointer to data container
 */
CScalarLoss::CScalarLoss(CModel* &model, CVecData* &data)
   : CLoss(model, 1, data->dim(), data->bias()),
     _data(data),
	 m(0),
	 f(0),
	 l(0)
{
   m = _data->slice_size();
   scalingFactor = 1.0/_data->size();
   f = new TheMatrix(m, 1, SML::DENSE);
   l = new TheMatrix(m, 1, SML::DENSE);  // it will be slower if we use sparse l	
}


/**  Destructor
 */
CScalarLoss::~CScalarLoss()
{
   if(f) delete f;
   if(l) delete l;
}


/**  Compute loss function value and gradient of loss function w.r.t. w
 *   
 *   @param loss [write] loss function value
 *   @param grad [write] gradient of loss function w.r.t w
 */
void CScalarLoss::ComputeLossAndGradient(double& loss, TheMatrix& grad)
{
   TheMatrix &w = _model->GetW();
   loss = 0;
   grad.Zero();
   _data->ComputeF(w, *f);   // f := X*w
   LossAndGrad(loss, *f, *l);
   loss = loss*scalingFactor;
   _data->Grad(*l, grad);
   grad.Scale(scalingFactor);
}


/** Compute loss function value given f = X*w
 *   
 *   @param loss [write] loss function value
 */
void CScalarLoss::ComputeLoss(double& loss)
{
   TheMatrix &w = _model->GetW();
   loss = 0; 
   _data->ComputeF(w, *f);   // f := X*w
   Loss(loss, *f);
   loss = loss*scalingFactor; 
}



/**  Evaluate the performance of the model on the examples in _data
 *
 *   @param w [read] weight vector
 */

void CScalarLoss::Evaluate(CModel* model)
{
   // sanity check
   assert(_data);
   
   TheMatrix &w = model->GetW();  
   
   // f := X*w
   _data->ComputeF(w, *f);
   TheMatrix predict(*f);
   Transform(predict);
   
   // evaluate
   Perf(*f, predict);
   
   // write outputs to files
   int m = f->Length();
   double* f_array = f->Data();   // pointer for fast access
   double* p_array = predict.Data();   // pointer for fast access
   
   Configuration &config = Configuration::GetInstance();
   bool outputPrediction = false;

   if(config.IsSet("Prediction.outputFvalAndLabels"))
      outputPrediction = config.GetBool("Prediction.outputFvalAndLabels");
   
   if(outputPrediction)
   {
      bool success = false;
      string predictedLabelsFn = "predictedLabels";
      if(config.IsSet("Prediction.predictedLabelsFile"))
         predictedLabelsFn = config.GetString("Prediction.predictedLabelsFile");
      
      success = WriteFile<double*>(predictedLabelsFn, m, p_array);
      if(success) cout << "Predicted labels file written." << endl;
      
      string decisionFunctionValuesFn = "decisionFunctionValuesFile";
      if(config.IsSet("Prediction.decisionFunctionValuesFile"))
         decisionFunctionValuesFn = config.GetString("Prediction.decisionFunctionValuesFile");
      
      success = WriteFile<double*>(decisionFunctionValuesFn, m, f_array);
      if(success) cout << "Decision function values file written." << endl;
   }
}

#endif
