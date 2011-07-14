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
 * Created: (12/01/2009) 
 *
 * Last Updated:
 */

#ifndef _FBETALOSS_HPP_
#define _FBETALOSS_HPP_

#include <vector>

#include "common.hpp"
#include "sml.hpp"
#include "binaryclassificationloss.hpp"
#include "model.hpp"

/**  Class to represent binary F-Beta loss function:
 *
 *   Reference:
 *   Thorsten Joachims, "A Support Vector Method for Multivariate Performance Measures", ICML 2005.
 */
class CFBetaLoss : public CBinaryClassificationLoss 
{
   protected:
      /** Number of positive training labels
       */
      unsigned int m_pos;
      
      /** Number of negative training labels
       */
      unsigned int m_neg;                
      
      /** Temporary indices
       */
      int* orig_pidx;
      int* orig_nidx;
      int* pidx; 
      int* nidx;

      /** Partial sums of f for positive and negative examples
       */
      double *fpos;
      double *fneg;

      /** The beta in F-beta score
       */
      double beta;

      void Loss(double& loss, TheMatrix& f);
      void LossAndGrad(double& loss, TheMatrix& f, TheMatrix& l);
      double Delta(unsigned int Tpos, unsigned int Tneg)
      {
         double fbeta = (1.0+beta*beta)*Tpos / (Tpos + m_neg - Tneg + beta*beta*m_pos);
         return 100.0*(1-fbeta);
      }

   public:    
      
      CFBetaLoss(CModel* &model, CVecData* &data);		
      virtual ~CFBetaLoss();
};

#endif
