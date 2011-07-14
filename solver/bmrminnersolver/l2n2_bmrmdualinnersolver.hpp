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
 * Created      : 20/11/2007
 * Last Updated :
 */

#ifndef _L2N2_BMRMDUALINNERSOLVER_HPP_
#define _L2N2_BMRMDUALINNERSOLVER_HPP_

#include "model.hpp"
#include "sml.hpp"
#include "bmrminnersolver.hpp"


namespace L2N2_BMRMDualInnerSolver
{
   const double ZERO_EPS = 1e-20;
}


/** 
 *   when \Omega = 0.5|w|_2^2, the dual of the problem can be solved instead of the primal.
 */
class CL2N2_BMRMDualInnerSolver : public CBMRMInnerSolver 
{
   protected:
      
      /** Dimensionality of x vector before update with new gradient and removal of old gradients
       */
      int prevDim;
      
      /** The variables (i.e. Lagragian multipliers)
       */
      double *x;
      
      /** Linear part of the objective function
       */
      double *f;
      
      /** Hessian matrix of the objective function
       */
      double *Q;
      
      /** Constraint matrix
       */
      double *a;
      
      /** RHS of the (in)equality constraints
       */
      double *b;
      
      /** Lower bound vector for the variables
       */
      double *l;
      
      /** Upper bound vector for the variables
       */
      double *u;
      
      /** Tolerance for optimization error
       */
      double tol;
      
      /** Gradient set
       */
      std::vector<TheMatrix*> gradientSet;
      
      /** Offsets set
       */
      std::vector<double> offsetSet;
      
      /** Type of gradient to be store in gradientSet (e.g., SPARSE or DENSE)
       *  [default: DENSE]
       */
      int gradType;
      
      /** Time stamps indicating the last iteration number a gradient is active.
       */
      std::vector<int> activeTimeStamp;
      
      /** Time stamps indicating the iteration number when a gradient entered the set
       */
      std::vector<int> enterTimeStamp;
      
      /** Max. number of consecutive iterations a Lagrange multiplier may
       *  remain zero before removal
       *  [default: 10]
       */
      int gradIdleAge;
      
      /** Maximum number of gradients to keep in gradientSet.
       *  Once the number of gradients exceeds this, oldest gradients 
       *  will be "aggregated" into one.
       *  [default: 1000]
       */
      int maxGradSetSize;
      
      /** Whether to remove ALL idling gradients exceeded gradIdleAge or just
       *  the oldest one.
       *  [default: false]
       */
      bool removeAllIdleGrad;
      
      /** Scaling factor for QP
       *  [default: 1.0]
       */
      double QPScale;
      
      /** index for aggregated gradient
       */
      int aggGradIdx;
      
    
      /** Use Diagonal Quadratic Form regularizer?
       */
      bool useDiagQuadForm;
      
      /** Diagonal Positive Definite Matrix
       */
      TheMatrix *B;
      
      /** Inverse of Diagonal PD matrix B
       */
      TheMatrix *invB;
      
      
      /** Solve QP
       */
      virtual void SolveQP()=0;
      
      
      /** Update the solver with new gradient
       */
      virtual void Update(TheMatrix& a, double b);
      
      
      /** Remove ALL idle gradients and return a position for insertion of new gradient
       */
      virtual int RemoveAllIdleGradients();
      
      
      /** Remove the laziest idle gradient adn return a position for insertion of new gradient
       */
      virtual int RemoveLaziestIdleGradient();
      
      
      /** Aggregate gradients and return a position for insertion of new gradient
       */
      virtual int AggregateGradients();
      
      
      /** Get solution and lower bound of empirical risk
       */
      virtual void GetSolution(TheMatrix& w, double &objval);
      
      
#ifdef QPDEBUG
      /** Routine to check the correctness of the quadratic matrix of the objective function
       */
      void MatrixCorrectnessCheck();
#endif
      
   public:     
      
      /** Constructor
       */
      CL2N2_BMRMDualInnerSolver(double lambda);      
      
      
      /** Destructor
       */
      virtual ~CL2N2_BMRMDualInnerSolver();
      
      
      /** Solve the problem
       */
      virtual void Solve(TheMatrix& w, TheMatrix& a, double loss, double &objval);
      
      
      /** Compute the value of regularizer
       *
       *  @param w [read] weight vector
       *  @return regularizer value
       */
      virtual double ComputeRegularizerValue(TheMatrix &w);
      
      
      /** Reset the gradientSet, offsetSet, and mem for QP
       */
      virtual void Reset();
      
      
      /** With good QP tolerance annealing heuristic 
       *  the whole problem can be solved in much less number of iterations
       */
      virtual void SetTolerance(const double &theTolerance) {tol = theTolerance;}
      
      /** 
       */
      virtual void ChangeLambda(double newLambda);
};

#endif
