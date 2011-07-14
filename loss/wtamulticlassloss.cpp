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

#ifndef _WTAMULTICLASSLOSS_CPP_
#define _WTAMULTICLASSLOSS_CPP_

#include "common.hpp"
#include "wtamulticlassloss.hpp"
#include "configuration.hpp"
#include <sstream>


using namespace std;


/** Constructor
 *
 *  @param w    [r/w] external pointer to weight vector
 *  @param data [r/w] pointer to data container
 *  @param additiveLabelLoss [read] flag to determine whether to "add" label loss to loss function or not (default: true)
 */
CWTAMulticlassLoss::CWTAMulticlassLoss(CModel* &model, CVecData* &data, bool additiveLabelLoss)
   : _data(data),
     f(0),
     matW(0),
     matG(0),
     g(0),
     _additiveLabelLoss(additiveLabelLoss),
     _indloss(0)
{
   // initialize model
   if(! model->IsInitialized())
   {
      if(! data->HasLabel())
         throw CBMRMException("Error: dataset contains no labels; unable to determine number of classes!",
                              "CWTAMulticlassLoss::CWTAMulticlassLoss()");
      numOfClass = (int)data->MaxLabel();
      model->Initialize(numOfClass, _data->dim(), _data->bias());
   }
   else
   {
      numOfClass = model->GetNumOfW();
   }
   
   // keep a pointer to model
   _model = model;
   
   
   // check if labels are suitable for multiclass classification
   if(_data->HasLabel() && ((int)_data->MinLabel() < 1.0)) 
   {
      throw CBMRMException("Labels must be positive integers",
                           "CWTAMulticlassLoss::CWTAMulticlassLoss()");
   }
   
   // check if we can access the examples as individual vectors
   if(! _data->HasMatrixRowView()) 
   {
      if(verbosity)
         cout << "CWTAMulticlassLoss: creating explicit views to examples... " << endl;
      _data->CreateFeatureMatrixRowViews();
      _data->CreateLabelMatrixRowViews();		
	}
   
   
   Configuration &config = Configuration::GetInstance();
   
   // default margin scaling function
   Gamma = &CWTAMulticlassLoss::MarginRescaling;
   marginScalingFunctionType = "MARGIN_RESCALING";
   
   if(config.IsSet("WTA_MULTICLASS.marginScalingFunctionType"))
   {
      if(config.GetString("WTA_MULTICLASS.marginScalingFunctionType") == "SLACK_RESCALING")
      {
         Gamma = &CWTAMulticlassLoss::SlackRescaling;
         marginScalingFunctionType = "SLACK_RESCALING";
      }
   }
   
   // default label loss function
   Delta = &CWTAMulticlassLoss::ZeroOne;
   labelLossType = "ZERO_ONE";
   
   if(config.IsSet("WTA_MULTICLASS.labelLossType"))
   {
      if(config.GetString("WTA_MULTICLASS.labelLossType") == "SQUARED_DIFFERENCE")
      {
         Delta = &CWTAMulticlassLoss::SquaredDifference;
         labelLossType = "SQUARED_DIFFERENCE";
      }
      if(config.GetString("WTA_MULTICLASS.labelLossType") == "ZERO_ONE_HUNDRED")
      {
         Delta = &CWTAMulticlassLoss::ZeroOneHundred;
         labelLossType = "ZERO_ONE_HUNDRED";
      }
      
   }
   
   // determine the number of classes from dataset
   m = _data->slice_size();
   scalingFactor = scalingFactor/_data->size();
   matW = new TheMatrix(numOfClass, _data->dim(), SML::DENSE);
   matG = new TheMatrix(numOfClass, _data->dim(), SML::DENSE);
   f = new TheMatrix(1, numOfClass, SML::DENSE);
   
   //i.e. g[i] is the i-th row of matG
   g = (TheMatrix**)malloc(sizeof(TheMatrix*)*numOfClass);
   for(unsigned int i=0; i<numOfClass; i++)
      g[i] = matG->CreateMatrixRowView(i);
   
   // prepare spaces for individual loss
   _indloss = new double[m];
   memset(_indloss, 0, sizeof(double)*m);
   
   if(verbosity > 0)
   {
      cout << "In CWTAMulticlassLoss::CWTAMulticlassLoss()" << endl;
      cout << "Problem: " << numOfClass << "-class classification" << endl;
      cout << "Margin scaling function type: " << marginScalingFunctionType << endl;
      cout << "Label loss type: " << labelLossType << endl << endl;                
   }
}


/**  Destructor 
 */
CWTAMulticlassLoss::~CWTAMulticlassLoss()
{
   if(matW) delete matW;
   if(matG) delete matG;
   if(f) delete f;
   if(g) 
   { 
      for(unsigned int i=0; i<numOfClass; i++) 
         delete g[i];
      free(g); 
      g = 0;
   }
   if(_indloss) delete _indloss;
}



/**  Compute loss
 */
void CWTAMulticlassLoss::ComputeLoss(double& loss)
{
   TheMatrix &w = _model->GetW();
   double y_i = 1;
   double ybar_i = 1;
   double f_y_i = 0;
   double f_ybar_i = 0;
   double loss_i = 0;
   double max_loss_i = 0;
   double labelloss_i = 0;
   double marginscale = 0.0;
   
   // reset slots for individual loss
   memset(_indloss, 0, sizeof(double)*m);
   
   // assign the content of w (1D) to matW (matrix)
   matW->Assign(w);  // chteo: to be replace with TheMatrix::StackRow()
   loss = 0;
   const TheMatrix &Y = _data->labels();
   
   for(unsigned int i=0; i<m; i++) 
   {
      _data->ComputeFi(*matW, *f, i);
      Y.Get(i, y_i);
      f->Get((int)y_i-1, f_y_i);
      max_loss_i = 0.0;
      loss_i = 0.0;
      labelloss_i = 0.0;
      marginscale = 0.0;
      ybar_i = -1;
      
      for(unsigned int j=0; j<numOfClass; j++) 
      {
         // ignore true label
         if((unsigned int)y_i != j+1)
         {
            f->Get(j, f_ybar_i);
            labelloss_i = (this->*Delta)((int)y_i, j+1);
            marginscale = (this->*Gamma)(labelloss_i);
            loss_i = marginscale * (f_ybar_i - f_y_i); // + labelloss
			
            // this is modified for ramp loss!
            if(_additiveLabelLoss)
               loss_i += labelloss_i;
            
            if(loss_i > max_loss_i) 
            {
               max_loss_i = loss_i;
               ybar_i = j+1;
            }
         }
      }
      
      if(max_loss_i > 0.0)
      {
         loss += max_loss_i;
         _indloss[i] = max_loss_i * scalingFactor;
      }
   }
   loss *= scalingFactor;  
}



/**   Compute loss and gradient
 *    Note: class label starts from 1
 */
void CWTAMulticlassLoss::ComputeLossAndGradient(double& loss, TheMatrix& grad)
{
   TheMatrix &w = _model->GetW();
   double y_i = 1;
   double ybar_i = 1;
   double f_y_i = 0;
   double f_ybar_i = 0;
   double loss_i = 0;
   double max_loss_i = 0;
   double labelloss_i = 0;
   double marginscale = 0;
   double max_marginscale = 0;
   
   // reset slots for individual loss
   memset(_indloss, 0, sizeof(double)*m);
   
   // assign the content of w (1D) to matW (matrix)
   matW->Assign(w);  
   matG->Zero();
   loss = 0;
   const TheMatrix &Y = _data->labels();
   int errcnt = 0;
   
   for(unsigned int i=0; i<m; i++) 
   {
      _data->ComputeFi(*matW, *f, i);
      Y.Get(i, y_i);
      f->Get((int)y_i-1, f_y_i);
      max_loss_i = 0.0;
      loss_i = 0.0;
      labelloss_i = 0.0;
      marginscale = 0.0;
      max_marginscale = 0.0;
      ybar_i = -1;
      
      
      // Do: argmax_y Gamma(Delta(y,ybar))(\f_ybar - f_y) + Delta(y,ybar), where f_* := <w,*>
      for(unsigned int j=0; j < numOfClass; j++) 
      {
         if((unsigned int)y_i != j+1)
         {
            f->Get(j, f_ybar_i);
            labelloss_i = (this->*Delta)((int)y_i, j+1);
            marginscale = (this->*Gamma)(labelloss_i);                               		
            loss_i = marginscale * (f_ybar_i - f_y_i); // + labelloss;
            
            // this is modified for ramp loss
            if(_additiveLabelLoss)
               loss_i += labelloss_i;
			
            if(loss_i >= max_loss_i) 
            {
               max_loss_i = loss_i;
               max_marginscale = marginscale;
               ybar_i = j+1;
            }
         }
      }
      
      if(max_loss_i > 0.0)
      {
         loss += max_loss_i;
         _data->AddElement(*g[(int)y_i-1], i, -max_marginscale);
         _data->AddElement(*g[(int)ybar_i-1], i, max_marginscale);
         errcnt++;
         _indloss[i] = max_loss_i * scalingFactor;
      }
   }
   
   grad.Assign(*matG);
   grad.Scale(scalingFactor);
   loss *= scalingFactor;
}



/**  Evaluate the performance of the model on the examples in _data
 *
 *   @param model [read] model object containing a weight vector
 */

void CWTAMulticlassLoss::Evaluate(CModel* model)
{
   // sanity check
   assert(model);
   assert(_data);
   assert(m > 0);
   
   if(not _data->HasLabel())
   {
      throw CBMRMException("Data object does not contain labels","CWTAMulticlassLoss::Evaluate()");
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
   int errorCnt = 0;
   double y_i = 0;;
   const TheMatrix &Y = _data->labels();
   
   for(unsigned int i=0; i < m; i++)
   {
      Y.Get(i, y_i);
      if(ybar[i] != (int)y_i) 
         errorCnt++;
   }
   
   double loss = 0.0;
   ComputeLoss(loss);
   
   // performance
   cout << "\nPerformance:" << endl;
   cout << "1. Accuracy (%): " << 100.0*(m-errorCnt)/m << endl;
   cout << "2. Misclassification: " << errorCnt << endl;
   cout << "3. Loss: " << loss << endl;
   
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


void CWTAMulticlassLoss::DecisionAndPrediction(CModel* model, int *ybar, double *f_ybar)
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


void CWTAMulticlassLoss::GetIndividualLoss(double* &indloss, unsigned int &len)
{
   indloss = _indloss;
   len = m;
}

#endif
