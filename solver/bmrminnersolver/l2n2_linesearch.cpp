
#ifndef _L2N2_LINESEARCH_CPP_
#define _L2N2_LINESEARCH_CPP_

#include "l2n2_linesearch.hpp"

using namespace std;


double CL2N2_LineSearch::ComputeRegularizerValue(TheMatrix &w)
{
	double regval = w.Norm2();
	regval = 0.5*bmrmLambda*regval*regval;
	return regval;
}


void CL2N2_LineSearch::Solve(TheMatrix& w, TheMatrix& a, double loss, double &objval)
{
    iter++;
    
    double slope = 0.0;
    double gradientNorm2 = 0.0;
    double w_dot_a = 0.0;
    double K = 0.0;
    double eta = 0.0;
    
    slope = bmrmLambda*prevWNorm2Sq + loss - r;  // \gamma_t
    a.Norm2(gradientNorm2);
    w.Dot(a,w_dot_a);
    K = bmrmLambda*bmrmLambda*prevWNorm2Sq + 2.0*bmrmLambda*w_dot_a + gradientNorm2*gradientNorm2;
    K /= bmrmLambda;   
    eta = min(1.0, slope/K);

    // make sure sum(alpha) == 1
    if(iter == 1) eta = 1.0;

    r = (1.0-eta)*r + eta*(loss - w_dot_a);
    w.Scale(1.0-eta);
    w.ScaleAdd(-eta/bmrmLambda, a);
    w.Norm2(prevWNorm2Sq);
    prevWNorm2Sq = prevWNorm2Sq*prevWNorm2Sq;
    objval = -bmrmLambda*0.5*prevWNorm2Sq + r;
}

#endif
