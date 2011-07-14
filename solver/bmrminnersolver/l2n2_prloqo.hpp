
 /* Purpose      : solves quadratic programming problem for pattern recognition
 *                and regression estimation for support vectors
 *
 * Note         : adapted from pr_loqo
 */

#ifndef _L2N2_PRLOQO_HPP_
#define _L2N2_PRLOQO_HPP_

//#include "common.hpp"
#include "sml.hpp"
#include "model.hpp"
#include "l2n2_bmrmdualinnersolver.hpp"
#include "configuration.hpp"

/* 
 * solve the quadratic programming problem
 *
 * minimize    c' * x + 1/2 x' * H * x
 * subject to  b <= A*x <= b+r
 *             l <= x <= u
 *
 *  for a documentation see R. Vanderbei, LOQO: an Interior Point Code
 *                          for Quadratic Programming
 *
 *  NOTE: the role of r and b are swapped to maintain
 *        consistency across different solvers.
 */

/*
 * n   : number of primal variables
 * m   : number of constraints (typically 1)
 * h_x : dot product matrix (n.n)
 * a   : constraint matrix (n.m)
 * b   : constant term (m)
 * l   : lower bound (n)
 * u   : upper bound (m)
 *
 * primal : workspace for primal variables, has to be of size 3 n
 *
 *  x = primal;			n
 *  g = x + n;			n
 *  t = g + n;			n
 *
 * dual : workspace for dual variables, has to be of size m + 2 n
 *
 *  y = dual;			m
 *  z = y + m;			n
 *  s = z + n;			n
 *
 * verb       : verbosity level
 * sigfig_max : number of significant digits
 * counter_max: stopping criterion
 * restart    : 1 if restart desired
 *
 */

namespace PRLOQO {
   enum PRLOQO_PHASE { PREDICTOR=1, CORRECTOR=2 };
   enum PRLOQO_VERBOSITY { QUIET=0, STATUS=4, FLOOD=5 };
   enum PRLOQO_MESSAGE { OPTIMAL_SOLUTION=1, SUBOPTIMAL_SOLUTION=2, ITERATION_LIMIT=3,
                         PRIMAL_INFEASIBLE=4, DUAL_INFEASIBLE=5, PRIMAL_AND_DUAL_INFEASIBLE=6,
                         INCONSISTENT=7, PRIMAL_UNBOUNDED=8, DUAL_UNBOUNDED=9, TIME_LIMIT=10 };   
}


class CL2N2_prLOQO : public CL2N2_BMRMDualInnerSolver
{
private:
    double sigfig_max;
    double initBound;  // bound for variables initialization
    int maxIntPointIter;
    int gradSetSizeIncrement;
    int cur_mem_n;  // mem size for # of primal var
    int cur_mem_m;  // mem size for # of dual var
    double *r;
    double *b_local;  // b is already used in superclass as RHS of singly linearly inequality constraint
    
    // mem used in Solve()
    double *primal;    // workspace for primal var
    double *dual;      // workspace for dual var
    double *workspace;
    double *diag_h_x;
    double *h_x;
    double *h_y;
    double *c_x;
    double *c_y;
    double *h_dot_x;
    double *rho;
    double *nu;
    double *tau;
    double *alpha;
    double *beta;
    double *sigma;
    double *gamma_z;
    double *gamma_w;  
    double *gamma_s;  
    double *gamma_q;  
    double *hat_nu;
    double *hat_tau;
    double *hat_alpha;
    double *hat_beta;
    double *delta_x;
    double *delta_y;
    double *delta_s;
    double *delta_z;
    double *delta_g;
    double *delta_t;
    double *delta_w;
    double *delta_p;
    double *delta_q;
    double *delta_v;
    double *d;
    double *e;
    
    
    // chteo: hx and hy_aug are redundant as we are doing a quick hack
    double *chol_hx;
    double *hxinv_a;
    double *P;
    double hy_aug;
    
    void choldc(double a[], int n, double p[]);
    void cholsb(double a[], int n, double p[], double b[], double x[]);
    void chol_forward(double a[], int n, double p[], double b[], double x[]);
    void chol_backward(double a[], int n, double p[], double b[], double x[]);
    void matrix_vector(int n, double m[], double x[], double y[]);
    void solve_reduced(int n, int m, double h_x[], double h_y[], 
		       double a[], double x_x[], double x_y[],
		       double c_x[], double c_y[],
		       double workspace[], int step);
    
    void reallocPrimalMem(int n);
    void reallocDualMem(int m);
    
//     /** Dump QP into files.
//      *  H.txt -- hessian matrix
//      *  c.txt -- linear part of the objective
//      */
//     void PrintQP();
 

protected:
    /** Solve QP
     */
    void SolveQP();
    
public:
    CL2N2_prLOQO(double lambda);
    
    virtual ~CL2N2_prLOQO();
    
    /** Lower tolerance causes dense interior point solution
     *  So, prLOQO solves every QP to a fixed sigfig.
     */
    virtual void SetTolerance(const double &theTolerance) 
	{
	    Configuration &config = Configuration::GetInstance();
	    if(config.IsSet("L2N2_prLOQO.annealMaxSigfig") && config.GetBool("L2N2_prLOQO.annealMaxSigFig"))
	    {
            sigfig_max = std::min(12.0, -log10(theTolerance*1e-1));
            //std::cout << "tol = " << theTolerance << "  sigfig_max = " << sigfig_max << std::endl;
	    }
	}
    
};




#endif
