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

#ifndef _MULTILABELLOSS_CPP_
#define _MULTILABELLOSS_CPP_

#include "common.hpp"  // WriteFile()
#include "multilabelloss.hpp"
#include "configuration.hpp"
#include <sstream>

using namespace std;

/** Constructor
 *
 *  @param w    [r/w] external pointer to weight vector
 *  @param data [r/w] pointer to data container
 *  @param additiveLabelLoss [read] flag to determine whether to "add" label loss to loss function or not (default: true)
 */
CMultilabelLoss::CMultilabelLoss(CModel* &model, CMultilabelVecData* &data)
   : _data(data),
     f(0),
     matW(0),
     matG(0),
     g(0)          
{
   // initialize model
   if(! model->IsInitialized())
   {
      if(! data->HasLabel())
         throw CBMRMException("Error: dataset contains no labels; unable to determine number of classes!",
                              "CMultilabelLoss::CMultilabelLoss()");
      numOfClass = (int)data->MaxLabel();
      model->Initialize(numOfClass, _data->dim(), _data->bias());
   }
   else
   {
      numOfClass = model->GetNumOfW();
   }
   
   // keep a pointer to model
   _model = model;
   
   // check if labels are suitable for multilabel loss
   if(_data->HasLabel() && (_data->MinLabel() < 1.0)) 
      throw CBMRMException("Labels must be positive integers","CMultilabelLoss::CMultilabelLoss()");
   
   // check if we can access the examples as individual vectors
   if(! _data->HasMatrixRowView()) 
   {
      if(verbosity)
         cout << "CMultilabelLoss: creating explicit views to examples... " << endl;
      _data->CreateFeatureMatrixRowViews();
   }	
   
   
   m = _data->slice_size();  // determine the number of classes from dataset
   scalingFactor = scalingFactor/_data->size();
   matW = new TheMatrix(numOfClass, _data->dim(), SML::DENSE);
   matG = new TheMatrix(numOfClass, _data->dim(), SML::DENSE);
   f = new TheMatrix(1, numOfClass, SML::DENSE);
   Ybar.resize(numOfClass, true);
   
   //i.e. g[i] is the i-th row of matG
   g = (TheMatrix**)malloc(sizeof(TheMatrix*)*numOfClass);
   for(unsigned int i=0; i<numOfClass; i++)
      g[i] = matG->CreateMatrixRowView(i);
   
   if(verbosity > 0)
   {
      cout << "In CMultilabelLoss::CMultilabelLoss()" << endl;
      cout << "Problem: " << numOfClass << "-class multilabel classification" << endl;
   }
}


/**  Destructor 
 */
CMultilabelLoss::~CMultilabelLoss()
{
   if(matW) {delete matW; matW = 0;}
   if(matG) {delete matG; matG = 0;}
   if(f) {delete f; f = 0;}
   if(g) 
   { 
      for(unsigned int i=0; i<numOfClass; i++)
         delete g[i];            
      free(g); 
      g = 0;
   }
}


/**  Compute loss
 */
void CMultilabelLoss::ComputeLoss(double& loss)
{
   TheMatrix dummygrad(_model->GetW());
   ComputeLossAndGradient(loss, dummygrad);
}


/**   Compute loss and gradient
 *    Note: class label starts from 1
 */
void CMultilabelLoss::ComputeLossAndGradient(double& loss, TheMatrix& grad)
{
   TheMatrix &w = _model->GetW();
   unsigned int y_i = 0;        
   unsigned int ybar_i = 0;
   double f_y_i = 0;
   double f_ybar_i = 0;
   
   // assign the content of w (1D) to matW (matrix)
   matW->Assign(w);  
   matG->Zero();
   loss = 0;
   const vector<vector<double> > &Y = _data->labels();	
   int errcnt = 0;
   
   for(unsigned int i=0; i<m; i++) 
   {
      // set Ybar for i-th example
      for(unsigned int j=0; j < Y[i].size(); j++)
         Ybar[(int)Y[i][j]-1] = false;
      
      y_i = 0;
      f_y_i = SML::INFTY;
      ybar_i = 0;
      f_ybar_i = -SML::INFTY;
      
      _data->ComputeFi(*matW, *f, i);
      
      
      // look for min_{y} f(x,y) and max_{ybar} f(x,ybar), where y is in true label set and y' in the complement   
      double f_k = 0;
      //unsigned int idx = 0;
      
      for(unsigned int k = 0; k < numOfClass; k++)
      {
         if(Ybar[k] == false)
         {
            Ybar[k] = true;  // reset for next use
            continue;
         }
         
         f->Get(k, f_k);
         if(f_k >= f_ybar_i)
         {
            f_ybar_i = f_k;
            ybar_i = k+1;
         }
      }
      
      for(unsigned int k = 0; k < Y[i].size(); k++)
      {
         f->Get((int)Y[i][k]-1, f_k);
         if(f_k <= f_y_i)
         {
            f_y_i = f_k;
            y_i = (unsigned int)Y[i][k];
         }
      }
      
      // loss := max(0, 1 + max_{y'} f(x,y') - min_{y} f(x,y)) 
      double loss_i = max(0.0, 1.0 + f_ybar_i - f_y_i);
      
      if(loss_i > 0.0)
      {
         loss += loss_i;
         // gradient for example i
         _data->AddElement(*g[y_i-1], i, -1.0);                
         _data->AddElement(*g[ybar_i-1], i, 1.0);
      }
      // Error only if ybar gets more weight than yi
      if(f_ybar_i > f_y_i)
         errcnt++;
   }
   
   grad.Assign(*matG);
   grad.Scale(scalingFactor);	
   loss *= scalingFactor;        
}

 


/**  Evaluate the performance of the model on the examples in _data
 *
 *   @param model [read] model object containing a weight vector
 */

void CMultilabelLoss::Evaluate(CModel* model)
{
   // sanity check
   assert(model);
   assert(_data);
   assert(m > 0);
   
   if(! _data->HasLabel())
   {
      throw CBMRMException("Data object does not contain labels","CMultilabelLoss::Evaluate()");
   }
   
   // more defensive precaution
   assert(_data->HasLabel());
   assert(_data->dim() == model->GetDimOfW());
   
   if(numOfClass != model->GetNumOfW())
   {
      if(matW) 
      {
         delete matW;                        
         matW = new TheMatrix(model->GetNumOfW(), model->GetDimOfW(), SML::DENSE);
      }                        
   }
   
   numOfClass = min(numOfClass, model->GetNumOfW());                      
   assert(matW);
   
   int* ybar = new int[m];
   double* f_ybar = new double[m];
   
   DecisionAndPrediction(model, ybar, f_ybar);
   
   // evaluation	
   int oneErrorCnt = m;
   const vector<vector<double> > &Y = _data->labels();
   
   for(unsigned int i=0; i < m; i++)
   {
      for(unsigned int j=0; j<Y[i].size(); j++)
         if(ybar[i] == Y[i][j])
         {
            oneErrorCnt--;
            break;
         }
   }
   
   // performance
   cout << "\nPerformance:" << endl;
   cout << "1. Acc@1 (%): " << 100.0*(m-oneErrorCnt)/m << endl;
   
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
      
      success = WriteFile<int*>(predictedLabelsFn, m, ybar);
      if(success) cout << "Predicted labels file written." << endl;
      
      string decisionFunctionValuesFn = "decisionFunctionValuesFile";
      if(config.IsSet("Prediction.decisionFunctionValuesFile"))
         decisionFunctionValuesFn = config.GetString("Prediction.decisionFunctionValuesFile");
      
      success = WriteFile<double*>(decisionFunctionValuesFn, m, f_ybar);
      if(success) cout << "Decision function values file written." << endl;
   }
   
   // clean up
   if(ybar) delete [] ybar;
   if(f_ybar) delete [] f_ybar;
}


void CMultilabelLoss::DecisionAndPrediction(CModel* model, int *ybar, double *f_ybar)
{
   // sanity check	
   assert(ybar);
   assert(f_ybar);
   
   TheMatrix &w = model->GetW();        
   matW->Assign(w);        	
	
   for(unsigned int i=0; i < m; i++)
   {
      _data->ComputeFi(*matW, *f, i);
      double max_f_ybar_i = -SML::INFTY;
      
      for(unsigned int j = 0; j < numOfClass; j++)
      {
         f->Get(j, f_ybar[i]);
         if(f_ybar[i] > max_f_ybar_i)
         {
            max_f_ybar_i = f_ybar[i];
            ybar[i] = j+1;
         }
      } 		
   }                        
}

#endif
