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
 * Created      : 25/07/2007
 * Last Updated :
 */

/* 
 * Purpose      : L1-norm Linear Program solver using COIN-OR Clp 
 *
 */

#ifndef _L1N1_CLP_HPP_
#define _L1N1_CLP_HPP_

#include "sml.hpp"
#include "bmrminnersolver.hpp"
#include <vector>

// from COIN-OR Clp
#include "ClpSimplex.hpp"

/* 
 * solve the 1-norm linear program
 *
 * minimize    c' * x + |x|_1
 * subject to  r <= A*x <= r+b
 *             l <= x <= u
 *
 */
class CL1N1_Clp : public CBMRMInnerSolver 
{
   protected:
      /** CLP simplex solver
       */
      ClpSimplex *sim;
      
      /** Scaling factor for LP objective 	
       */
      double LPScale;
      
      /** Pre-allocated indices array for new constaint
       */
      int *newRowIndices;
      
      /** Pre-allocated values array for new constraint
       */
      double *newRowElements;
      
      /** Latest time stamp for a constraint being active (i.e. getRowStatus(i) == ClpSimpex::atUpperBound)
       */
      std::vector<int> activeTimeStamp;
      
      /** Max. number of consecutive iterations a constraint may
       *  remain inactive (i.e. getRowStatus(i) == ClpSimpex::basic)
       *  [default: 10]
       */
      int gradIdleAge;
      
      /** If true then removel all inactive constraint (or gradient),
       *  otherwise, remove the oldest one.
       *  [default: false]
       */
      bool removeAllIdleGrad;
      
      /** Update the constraint matrix
       */
      virtual void Update(TheMatrix& a, double b);
      
      /** Return the updated solution
       */
      virtual void GetSolution(TheMatrix& w, double &objval);
      
   public:      
      
      /** Constructor
       */
      CL1N1_Clp(double lambda, const int &thedim);
      
      /** Destructor
       */
      virtual ~CL1N1_Clp();
      
      /** Solve the problem
       */
      virtual void Solve(TheMatrix& w, TheMatrix& a, double loss, double &objval);
      
      /** Compute the value of the regularizer
       *
       *  @param w [read] weight vector
       *  @return the regularizer value
       */
      virtual double ComputeRegularizerValue(TheMatrix &w);
      
      
      /** Clear constraint matrix
       */
      virtual void Reset();
      
      /** Change regularization constant lambda
       */
      virtual void ChangeLambda(double newLambda);
};

#endif
