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
 * Created: 16/05/2008
 *
 * Last Updated: 16/01/2009
 */

#include "linesearch_hingeloss.hpp"
#include "configuration.hpp"

using namespace std;


CLinesearch_HingeLoss::CLinesearch_HingeLoss(CModel* &model, CVecData* &data)
   : CBinaryClassificationLoss(model, data),
     ls(0),
     wb(0),
     fb(0),
     df(0),
     bmrmLambda(1.0),
     theta(0.1),
     loss_of_wb(0.0),
     eta(0.0)
{	
   Configuration &config = Configuration::GetInstance();
   bmrmLambda = config.GetDouble("BMRM.lambda");
   if(config.IsSet("Linesearch.theta"))
      theta = config.GetDouble("Linesearch.theta");
   
   unsigned int row = 0, col = 0;
   _model->GetW().Shape(row,col);
   
   wb = new TheMatrix(_model->GetW());
   p  = new TheMatrix(row,col,SML::DENSE);
   fb = new TheMatrix(1,m,SML::DENSE);
   df = new TheMatrix(1,m,SML::DENSE);
   delta = new TheMatrix(1,m,SML::DENSE);
   ls = new CHinge_Linesearch(model,this,data);

   // needed for AddElements
   _data->CreateFeatureMatrixRowViews();
}


CLinesearch_HingeLoss::~CLinesearch_HingeLoss() 
{
   if(wb) delete wb;
   if(p)  delete p;
   if(fb) delete fb;
   if(df) delete df;
   if(delta) delete delta;
}


void CLinesearch_HingeLoss::ComputeLoss(double &loss)
{
   TheMatrix grad(_model->GetW());
   ComputeLossAndGradient(loss, grad);
}


void CLinesearch_HingeLoss::ComputeLossAndGradient(double &loss, TheMatrix &grad)
{
   // todo:
   // 1. O(m) time to check if p is a descent direction. if so, do line-search, otherwise, skip it.
   
   eta = 0.0;
   TheMatrix &w = _model->GetW();
   
   // Y.*<w_t,X>
   _data->ComputeF(w, *f);
   f->ElementWiseMult(_data->labels());
   
   // Y.*<w_t-w_t^b,X>
   df->Assign(*f);
   df->Minus(*fb);
   
   // possible descent direction
   p->Assign(w);
   p->Minus(*wb);
   
   // line search
   //eta = ls->Linesearch(bmrmLambda, *wb, *p, scalingFactor, *fb, *df);
   eta = ls->Linesearch(*p, *fb, *df, bmrmLambda, *delta);
   
   // compute w_{t+1}^b and Y.*<w_{t+1}^b,X>
   if(eta > 0.0)
   {
      wb->ScaleAdd(eta, *p);
      fb->ScaleAdd(eta, *df);
   }
   
   
   // update w in model object which is actually used to compute loss and gradient
   w.Scale(theta);
   w.ScaleAdd((1-theta), *wb);
   
   // Y.*<w_{t+1}^c,X>
   f->Scale(theta);
   f->ScaleAdd((1-theta), *fb);  
   
   // compute loss and gradient using w_{t+1}^c
   loss = 0.0;
   grad.Zero();
   l->Zero();
   double* f_array = f->Data();  // pointer to memory location of f (faster element access)
   unsigned int len = f->Length();
   const TheMatrix &Y = _data->labels();
   double* Y_array = Y.Data();
   
   for(unsigned int i=0; i < len; i++) 
   {
      if(f_array[i] < 1.0) 
      {
         loss += 1-f_array[i];
         l->Set(i,-Y_array[i]);
      }
   }
   _data->Grad(*l,grad);
   loss *= scalingFactor;
   grad.Scale(scalingFactor);
   
   
   // compute loss using w_{t+1}^b
   int errcnt = 0;
   loss_of_wb = 0.0;
   f_array = fb->Data();
   for(unsigned int i=0; i < len; i++) 
   {
      if(f_array[i] < 1.0) 
      {
         loss_of_wb += 1-f_array[i];
         errcnt++;
      }
   }
   loss_of_wb *= scalingFactor;	
   //printf("errorrate=%.8e  eta=%.8e  ",100.0*errcnt/len, eta);	
}
