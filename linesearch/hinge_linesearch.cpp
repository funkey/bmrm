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
 * Authors: Jin Yu (jin.yu@nicta.com.au) 
 *          Simon Guenter
 *          Choon Hui Teo (choonhui.teo@anu.edu.au)
 *
 * Created: (14/02/2008)
 *
 * Last Updated: 16/01/2009
 */

#include "hinge_linesearch.hpp"

using namespace std;

CHinge_Linesearch::CHinge_Linesearch(CModel* model, CLoss* loss, CData *data) 
   : CLinesearch(model, loss, data) 
{
   alpha = vector<ScalarWithIndex> (numElements);
}

double CHinge_Linesearch::Linesearch(const TheMatrix &p, 
                                     TheMatrix &f,
                                     TheMatrix &df,
                                     const double &lambda, 
                                     TheMatrix &delta) 
{
   hinges.clear();
   
   const TheMatrix &w = _model->GetW();
   
   double b = 0.0;
   double h = 0.0;
   
   p.Dot(w, b);
   b *= lambda;
   
   p.Norm2(h);
   h = h * h;
   h *= lambda;
   
   // find hinges
   int hpi = 0;
   double min_hp = SML::INFTY, max_hp = -SML::INFTY;
   for (int i = 0; i < numElements; i++) 
   {
      double val1;
      double val2;
      f.Get(i, val1);
      df.Get(i, val2);
      if (val2 != 0) {
         val2 = (1 - val1) / val2;
         if (val2 > 0) {
            alpha[hpi].index = i;
            alpha[hpi].value = val2;
            if (val2 < min_hp) {
               min_hp = val2;
            }
            if (val2 > max_hp) {
               max_hp = val2;
            }
            
            hpi++;
         }
         
      }
   }
   
   // first guess of upper bound of eta
   double max_stepsize = max(2.0, 2 * min_hp);
   int firstposindex = 0;
   int lastposindex = hpi - 1;
   bool first = true;
   // supremum of subgradient
   double sup_subgrad = -1;
   double rho = 0.0, prev_rho = 0.0, hval = 0.0;
   eta = 0.0;
   
   while (hpi > 0) 
   {
      max_stepsize = min(max_stepsize, max_hp);
      
      for (; lastposindex >= firstposindex; lastposindex--) 
      {
         double thisval = alpha[lastposindex].value;
         if (thisval < max_stepsize)
            break;
      }
      
      for (int i = lastposindex; i >= firstposindex; i--) {
         double thisval = alpha[i].value;
         if (thisval > max_stepsize) {
            alpha[i].swap(alpha[lastposindex]);
            lastposindex--;
         }
      }

      //sort hinges in non-descending order
      sort(alpha.begin() + firstposindex, alpha.begin() + lastposindex + 1);
      
      if (first) {
         first = false;
         TheMatrix temp(f);
         //some value in-between 0 and the first hinge
         hval = alpha[0].value / 2.0;
         temp.ScaleAdd(hval, df);
         //initialize the indicator vector, delta
         for (int j = 0; j < numElements; j++) {
            double val;
            temp.Get(j, val);
            val = (val < 1);
            delta.Set(j, val);
         }
         
         //calculate <delta, \Delta f>
         delta.Dot(df, rho);
         // (sub)gradient at 0
         rho = C*rho - b;
         sup_subgrad = -rho;
      }

      int i = firstposindex;
      
      // find a hinge at which the supremum of the subgradient is non-negative
      while (sup_subgrad < 0) {
         prev_rho = rho;
         if (i == lastposindex + 1) {
            eta = SML::INFTY;
            break;
         } else {
            eta = alpha[i].value;
         }
         //update <delta, \Delta f>
         hval = alpha[i].value;
         while (alpha[i].value == hval && i<=lastposindex) {
            int index = alpha[i].index;
            double val;
            double dfval;
            delta.Get(index, val);
            df.Get(index, dfval);
            //update indicator vector delta and rho
            if (val == 0.0) {
               rho = rho + C*dfval;
               delta.Set(index, 1.0);
            } else {
               rho = rho - C*dfval;
               delta.Set(index, 0.0);
            }
            
            i++;
         }
         
         //calculate supremum of subdifferential at ith hinge
         sup_subgrad = hval * h - rho;
      }
      
      if (sup_subgrad >= 0 || lastposindex == hpi - 1) {
         // end because the supremum of subgradient
         // is positive or all hinges are tested.
         double proposedEta = prev_rho / h;
         // exact step size is min(eta, prev_rho / h);
         if (eta > proposedEta) {
            eta = proposedEta;
            //printf("exact step: %e\n", eta);
         }
         
         break;
      } else {
         // optimal eta not found, increase max step size
         max_stepsize = 2.0 * max_stepsize;
         firstposindex = lastposindex + 1;
         lastposindex = hpi - 1;
      }
	}

	return eta;
}
