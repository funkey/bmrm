
#ifndef _L2N2_LINESEARCH_HPP_
#define _L2N2_LINESEARCH_HPP_

#include "model.hpp"
#include "sml.hpp"
#include "bmrminnersolver.hpp"


/** 
 *   Solve the bmrm problem in the "primal" using linear search (as in the convergence proof)
 *   when \Omega(w) = 0.5*|w|_2^2
 */
class CL2N2_LineSearch : public CBMRMInnerSolver 
{
   protected:
      
      /** cache storing the value of b^T\alpha
       */
      double r;

      /** cache storing |w|_2^2 at previous iteration
       */
      double prevWNorm2Sq;
      
   public:     
      
      /** Constructor
       */
      CL2N2_LineSearch(double lambda) : CBMRMInnerSolver(lambda), r(0.0), prevWNorm2Sq(0.0) 
      {
	 if(verbosity > 0)
	    std::cout << "CL2N2_LineSearch instantiated!" << std::endl;
      }      

      /** Destructor
       */
      virtual ~CL2N2_LineSearch() {}
      
      /** Solve the problem
       */
      virtual void Solve(TheMatrix& w, TheMatrix& a, double loss, double &objval);

      virtual double ComputeRegularizerValue(TheMatrix &w);
      
      virtual void ChangeLambda(double newLambda) {bmrmLambda = newLambda;}
};

#endif
