/* Copyright (c) 2009 NICTA
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
 * Created: (14/11/2007) 
 *
 * Last Updated: (07/01/2009)   
 */

#ifndef _BMRM_HPP_
#define _BMRM_HPP_

#include <string>
#include "common.hpp"
#include "solver.hpp"
#include "bmrminnersolver.hpp"
#include "model.hpp"
#include "loss.hpp"
#include "sml.hpp"
#include "timer.hpp"

/**   Class for BMRM solver.
 *    This type of solver iteratively builds up a convex lower-bound of the 
 *      objective function, and performs minimization on the lower-bound.
 */
class CBMRM : public CSolver
{
   public:
      /** Constructor
       *  @param model [r/w] model object
       *  @param loss [read] loss object
       *  @param _procID [read] process ID (for parallelized bmrm)
       *  @param _nProc [read] number of processes (for parallelized bmrm)
       */
      CBMRM(CModel* model, CLoss* loss, int _procID=0, int _nProc=1);
      
      // Destructor
      virtual ~CBMRM();
      
      // Methods
      virtual void Train();
      
   protected:
      /** Verbosity level
       */
      int verbosity;
      
      /** Maximum number of BMRM iteration
       */
      unsigned int maxNumOfIter;
      
      /** Tolerance for epsilon termination criterion.
       *  epsilon = min_{t' <= t} J(w_{t'}) - J_t(w_t),
       *  where J_t is piece-wise linear approx. of J
       */
      double epsilonTol;
      
      /** Relative epsilonTol.
       *  Relative epsilon = epsilon / min_{t' <= t} J(w_{t'})
       */
      double relEpsilonTol;
      
      /**  Tolerance for gamma termination criterion.
       *   gamma = J(w_t) - J_t(w_t)
       */
      double gammaTol;
      
      /** Relative gammaTol.
       *  relGamma = gamma / J(w_t)
       */
      double relGammaTol;
      
      /** Regularization constant
       */
      double lambda;
      
      /** Prefix for intermediate/checkpoint model files
       */
      std::string checkpointPrefix;
      
      /** The number of iterations before saving a checkpoint model
       */
      unsigned int checkpointInterval;
      
      /** Selected type of checkpoint.
       */
      unsigned int checkpointMode;   
      
      /** Types of checkpoint mode.
       *  KEEP_ALL -- Keep a checkpoint model after every $checkpointInterval$ 
       *              iterations
       *  KEEP_LATEST -- Keep only the latest model
       */
      enum CHECKPOINT_MODE {KEEP_ALL, KEEP_LATEST};
      
      /* Pointer to inner solver e.g. qp solver
       */
      CBMRMInnerSolver* innerSolver;

      /** Process ID and number of processes for use in parallelized bmrm
       */
      int procID;
      int nProc;
      
      
      /** Validate all provided program parameters are good to use.
       *  E.g. lambda must be > 0.
       */
      virtual void ConfirmProgramParameters();

      
      /** Display iteration information
       */
      virtual void DisplayIterationInfo(unsigned int iter, double exactObjVal, 
                                        double approxObjVal, double epsilon, 
                                        double gamma, double loss, 
                                        double regVal, double curTime);
      

      /** Save model at every #checkpointInterval iterations
       */
      virtual void SaveCheckpointModel(unsigned int iter);
         
      
      /** Termination criteria check
       *
       */
      virtual int CheckTermination(unsigned int iter, double epsilon, double gamma, 
                                   double minExactObjVal, double exactObjVal);

      
      /** Adjust optimization tolerance of inner solver
       */
      virtual void AdjustInnerSolverOptTol(double& innerSolverTol, double prevEpsilon,
                                           double epsilon);

      /** Display information after training is done
       */
      virtual void  DisplayAfterTrainingInfo(unsigned int iter, double finalExactObjVal, 
                                             double approxObjVal, double loss, 
                                             TheMatrix& w_best, CTimer& lossAndGradientTime,
                                             CTimer& innerSolverTime, CTimer& totalTime);
      
};

#endif
