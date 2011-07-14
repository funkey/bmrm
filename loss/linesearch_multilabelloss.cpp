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
 *          Jin Yu (jin.yu@nicta.com.au)
 *
 * Created: 04/06/2008
 *
 * Last Updated: 11/10/2008
 */

#ifndef _LINESEARCH_MULTILABELLOSS_CPP_
#define _LINESEARCH_MULTILABELLOSS_CPP_

#include "linesearch_multilabelloss.hpp"
#include "multilabel_linesearch.hpp"
#include "configuration.hpp"

#define PROBLEM_SIZE_LIMIT 1e9

using namespace std;


CLinesearch_MultilabelLoss::CLinesearch_MultilabelLoss(CModel* &model,
                                                       CMultilabelVecData* &data) 
   : CMultilabelLoss(model, data), 
     F(0), 
     F_i(0), 
     Fb(0), 
     DF(0), 
     Wb(0), 
     P(0),
     ls(0), 
     bmrmLambda(1.0), 
     theta(0.1), 
     loss_of_wb(0.0), 
     eta(0.0) 
{
   if (double(numOfClass) * m >= PROBLEM_SIZE_LIMIT) 
   {
      throw CBMRMException(
         "Problem is too large; it will take large mem to keep 3e9 doubles!",
         "CLinesearch_MultilabelLoss()");
   }
   
   F = new TheMatrix(m, numOfClass, SML::DENSE);
   F_i = new TheMatrix();
   Fb = new TheMatrix(m, numOfClass, SML::DENSE);
   DF = new TheMatrix(m, numOfClass, SML::DENSE);
   Wb = new TheMatrix(1, numOfClass * _data->dim(), SML::DENSE);
   P = new TheMatrix(1, numOfClass * _data->dim(), SML::DENSE);
   ls = new CMultilabel_Linesearch(model, this, data);
   
   Configuration &config = Configuration::GetInstance();
   bmrmLambda = config.GetDouble("BMRM.lambda");
   if (config.IsSet("Linesearch.theta")) {
      theta = config.GetDouble("Linesearch.theta");
   } else {
      theta = 0.1;
   }
   
   if (config.IsSet("Linesearch.method")) {
      ls_method = config.GetInt("Linesearch.method");
   } else {
      ls_method = EXACT_LS;
   }
   
   // initialization of Wb and Fb
   matW->Assign(_model->GetW());
   for (unsigned int i = 0; i < m; i++) {
      Fb->CreateMatrixRowView(i, F_i);
      _data->ComputeFi(*matW, *F_i, i);
   }
   Wb->Assign(*matW);
}


CLinesearch_MultilabelLoss::~CLinesearch_MultilabelLoss() 
{
   if (F) delete F;
   if (F_i)	delete F_i;
   if (Fb) delete Fb;
   if (DF) delete DF;
   if (Wb) delete Wb;
   if (P) delete P;
   if (ls) delete ls;
}


/**   Compute loss and gradient
 *    Note: class label starts from 1
 */
void CLinesearch_MultilabelLoss::ComputeLossAndGradient(double& loss,
                                                        TheMatrix& grad) 
{
   eta = 0.0;
   TheMatrix &w = _model->GetW();
   bool calculate_loss = true, check_descent = true;
   
   // assign the content of w (1D) to matW (matrix)
	matW->Assign(w);
    
	// Compute <w,X>
	for (unsigned int i = 0; i < m; i++) 
    {
       F->CreateMatrixRowView(i, F_i);
       _data->ComputeFi(*matW, *F_i, i);
	}
    
	// Compute <w-wb,X> = <w,X> - <wb,X>
	DF->Assign(*F);
	DF->Minus(*Fb);

	// Compute possible descent direction, w-wb
	P->Assign(w);
	P->Minus(*Wb);

	// do line-search, and ask line search to calculate loss
	switch (ls_method) {
       case EXACT_LS:
          eta = ls->Linesearch(*Wb, *P, *Fb, *DF, bmrmLambda, calculate_loss,
                               check_descent);
          break;
       case INEXACT_LS:
          eta = ls->InexactLinesearch(*Wb, *P, *Fb, *DF, bmrmLambda, 1.0, 1e-3);
          break;
	default:
       printf("line search option unknown! please choose 0 for exact ls, 1 for inexact ls using cached slopes and offsets of active linear lines.n");
       exit(0);
    }
    
	// compute W_{t+1}^b and <W_{t+1}^b,X>
	double cur_theta = theta;
	if (eta > 0.0) {
       Wb->ScaleAdd(eta, *P);
       Fb->ScaleAdd(eta, *DF);
	}
    
	// update W as W_{t+1}^c
	w.Scale(theta);
	w.ScaleAdd((1 - cur_theta), *Wb);
    
	// compute <w_{t+1}^c,X>
	F->Scale(cur_theta);
	F->ScaleAdd((1 - cur_theta), *Fb);
    
	// compute loss and gradient at w_{t+1}^c
	const vector<vector<double> > &Y = _data->labels();
	loss = 0;
	matG->Zero();
	grad.Zero();
    
	unsigned int y_i = 0;
	unsigned int ybar_i = 0;
	double f_y_i = 0;
	double f_ybar_i = 0;
	double loss_i = 0;

	for (unsigned int i = 0; i < m; i++) {
       // set Ybar for i-th example
       for (unsigned int j = 0; j < Y[i].size(); j++)
          Ybar[(int) Y[i][j] - 1] = false;
       
       loss_i = 0;
       y_i = 0;
       f_y_i = SML::INFTY;
       ybar_i = 0;
       f_ybar_i = -SML::INFTY;
       
       F->CreateMatrixRowView(i, F_i);
       
       // look for min_{y} f(x,y) and max_{ybar} f(x,ybar), where y is in
       // true label set and ybar in the complement
       double f_k = 0;
       
       for (unsigned int k = 0; k < numOfClass; k++) {
          if (Ybar[k] == false) {
             Ybar[k] = true; // reset for next use
             continue;
          }
          
          F_i->Get(k, f_k);
          if (f_k > f_ybar_i) {
             f_ybar_i = f_k;
             ybar_i = k + 1;
          }
       }
       
       for (unsigned int k = 0; k < Y[i].size(); k++) {
          F_i->Get((int) (Y[i][k]) -1, f_k);
          if (f_k < f_y_i) {
             f_y_i = f_k;
             y_i = (unsigned int) Y[i][k];
          }
       }
       
       loss_i = max(0.0, 1.0 + f_ybar_i - f_y_i);
       
       if (loss_i > 0.0) {
          loss += loss_i;
          _data->AddElement(*g[y_i - 1], i, -1.0);
          _data->AddElement(*g[ybar_i - 1], i, 1.0);
       }
	}
    
	loss *= scalingFactor;
	grad.Assign(*matG);
	grad.Scale(scalingFactor);
    
	// compute loss using wb
    
	loss_of_wb = ls->GetLoss();
	//printf("eta=%.4e  ", eta);
}


/**  Compute loss at the given iterate w
 */
double CLinesearch_MultilabelLoss::ComputeLoss(TheMatrix &w) 
{
   unsigned int y_i = 0;
   unsigned int ybar_i = 0;
   double f_y_i = 0;
   double f_ybar_i = 0;
   
   // assign the content of w (1D) to matW (matrix)
   matW->Assign(w);
   double loss = 0;
   const vector<vector<double> > &Y = _data->labels();
   
   int errcnt = 0;
   
   for (unsigned int i = 0; i < m; i++) 
   {
      // set Ybar for i-th example
      for (unsigned int j = 0; j < Y[i].size(); j++)
         Ybar[(int) Y[i][j] - 1] = false;
      
      y_i = 0;
      f_y_i = SML::INFTY;
      ybar_i = 0;
      f_ybar_i = -SML::INFTY;
      
      _data->ComputeFi(*matW, *f, i);
      
      // look for min_{y} f(x,y) and max_{ybar} f(x,ybar), where y is in true label set and y' in the complement
      double f_k = 0;
      //unsigned int idx = 0;
      
      for (unsigned int k = 0; k < numOfClass; k++) {
         if (Ybar[k] == false) {
            Ybar[k] = true; // reset for next use
            continue;
         }
         
         f->Get(k, f_k);
         if (f_k >= f_ybar_i) {
            f_ybar_i = f_k;
				ybar_i = k + 1;
			}
		}

		for (unsigned int k = 0; k < Y[i].size(); k++) {
			f->Get((int) Y[i][k] - 1, f_k);
			if (f_k <= f_y_i) {
				f_y_i = f_k;
				y_i = (unsigned int) Y[i][k];
			}
		}

		double loss_i = max(0.0, 1.0 + f_ybar_i - f_y_i);

		if (loss_i > 0.0) {
			loss += loss_i;
		}
		// Error only if ybar gets more weight than yi
		if (f_ybar_i > f_y_i)
			errcnt++;
	}

	loss *= scalingFactor;

	return loss;
}


#endif
