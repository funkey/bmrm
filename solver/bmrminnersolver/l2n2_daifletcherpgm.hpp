/* Copyright (c) 2009, NICTA
 * All rights reserved.
 *
 * The contents of this file are subject to the Mozilla Public License Version
 * 1.1 (the "License"); you may not use this file except in compliance with
 * the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS IS" basis,
 * WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License
 * for the specific language governing rights and limitations under the
 * License.
 *
 * Authors      : Choon Hui Teo (ChoonHui.Teo@anu.edu.au)
 * Created      : 16/11/2006
 * Last Updated : 25/07/2007
 */

/* Notes :
 *
 * Dai-Fletcher Projected Gradient method for SVM [1].
 * Modified from the GPDT software [2] for inequality constraints.
 *
 * References:
 *
 *   [1] Y. H. Dai and R. Fletcher,
 *       New algorithms for singly linearly constrained quadratic programs
 *       subject to lower and upper bounds, 
 *       Math. Program., 2006.
 *
 *   [2] L. Zanni, T. Serafini, and G. Zanghirati,
 *       Parallel Software for Training Large Scale Support Vector Machines
 *       on Multiprocessor Systems,
 *       JMLR 7, 2006.
 */

#ifndef _L2N2_DAIFLETCHERPGM_HPP_
#define _L2N2_DAIFLETCHERPGM_HPP_


#include "sml.hpp"
#include "model.hpp"
#include "l2n2_bmrmdualinnersolver.hpp"
#include <iostream>
#include <cassert>


namespace DaiFletcherPGM {
   const double alpha_min = 1e-10;
   const double alpha_max = 1e10;
   const double EPS_SV = 1e-15;
   const double eps = 1e-20;
   const double tol_lam = 1e-15;
   const double tol_r = 1e-16;
}

class CL2N2_DaiFletcherPGM : public CL2N2_BMRMDualInnerSolver 
{
   public:      
      CL2N2_DaiFletcherPGM(double lambda);      
      virtual ~CL2N2_DaiFletcherPGM();
      
      /** Solve the QP
       */
      virtual void SolveQP();
      
   private: 
      int maxProjIter;
      int maxPGMIter;
      int *ipt, *ipt2, *uv;
      double *g, *y, *tempv, *d, *Qd, *t, *xplus, *tplus, *sk, *yk;
      int *flag;
      
      int ProjectDF(int n, double *a, double b, double *c, double *l, 
                    double *u, double *x, double &lam_ext);
      
      double ProjectR(double *x, int n, double lambda, double *a, 
                      double b, double *c, double *l, double *u);
                  
};

#endif
