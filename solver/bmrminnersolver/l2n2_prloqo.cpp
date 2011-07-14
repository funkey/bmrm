
 /* 
 * Purpose      : solves quadratic programming problem
 *
 * Note         : adapted from pr_loqo
 */

#ifndef _L2N2_PRLOQO_CPP
#define _L2N2_PRLOQO_CPP

#ifdef USE_LAPACK_ATLAS
extern "C" {
#include <clapack.h>
#include <cblas.h>
}
#endif

#include <fstream>
#include <cmath>
#include <cassert>
#include "l2n2_prloqo.hpp"
#include "configuration.hpp"

#define LAPACK_ERR_CHK(m,f,i) if((i)!=0) {cout<<(m)<<": "<<(f)<<" failed. info="<<(i)<<endl;exit(EXIT_FAILURE);}

using namespace std;

CL2N2_prLOQO::CL2N2_prLOQO(double lambda)
   : CL2N2_BMRMDualInnerSolver(lambda),
     sigfig_max(7.0),
     initBound(10.0),
     maxIntPointIter(5000),
     gradSetSizeIncrement(50),
     cur_mem_n(1),
     cur_mem_m(1),
     r(0),
     b_local(0),
     primal(0),
     dual(0),
     workspace(0),
     diag_h_x(0),
     h_x(0),
     h_y(0),
     c_x(0),
     c_y(0),
     h_dot_x(0),
     rho(0),
     nu(0),
     tau(0),
     alpha(0),
     beta(0),
     sigma(0),
     gamma_z(0),
     gamma_w(0),
     gamma_s(0),
     gamma_q(0),
     hat_nu(0),
     hat_tau(0),
     hat_alpha(0),
     hat_beta(0),
     delta_x(0),
     delta_y(0),
     delta_s(0),
     delta_z(0),
     delta_g(0),
     delta_t(0),
     delta_w(0),
     delta_p(0),
     delta_q(0),
     delta_v(0),
     d(0),
     e(0),
     chol_hx(0),
     hxinv_a(0),
     P(0)
{
    // synchronize members with configuration file
    Configuration &config = Configuration::GetInstance();
    
    if(config.IsSet("L2N2_prLOQO.maxIntPointIter"))
	maxIntPointIter = config.GetInt("L2N2_prLOQO.maxIntPointIter");
    
    if(config.IsSet("L2N2_prLOQO.maxSigFig"))  
	sigfig_max = config.GetDouble("L2N2_prLOQO.maxSigFig");
    
    // pre-allocate mem needed by solver
    reallocPrimalMem(cur_mem_n);
    
    // chote: remove this later as we assumed that there is only 1 linear constraint
    reallocDualMem(cur_mem_m);

    // scale QP problem by a factor of bmrmLambda
    //QPScale = bmrmLambda;
}



CL2N2_prLOQO::~CL2N2_prLOQO()
{
   // free memory
  if(r)         free(r);
  if(b_local)   free(b_local);
  if(primal)    free(primal);
  if(dual)      free(dual);
  if(workspace) free(workspace);
  if(diag_h_x)  free(diag_h_x);
  if(h_x)       free(h_x);
  if(h_y)       free(h_y);
  if(c_x)       free(c_x);
  if(c_y)       free(c_y);
  if(h_dot_x)   free(h_dot_x);
  if(rho)       free(rho);
  if(nu)        free(nu);
  if(tau)       free(tau);
  if(alpha)     free(alpha);
  if(beta)      free(beta);
  if(sigma)     free(sigma);
  if(gamma_z)   free(gamma_z);
  if(gamma_w)   free(gamma_w);
  if(gamma_s)   free(gamma_s);
  if(gamma_q)   free(gamma_q);
  if(hat_nu)    free(hat_nu);
  if(hat_tau)   free(hat_tau);
  if(hat_alpha) free(hat_alpha);
  if(hat_beta)  free(hat_beta);
  if(delta_x)   free(delta_x);
  if(delta_y)   free(delta_y);
  if(delta_s)   free(delta_s);
  if(delta_z)   free(delta_z);
  if(delta_g)   free(delta_g);
  if(delta_t)   free(delta_t);
  if(delta_w)   free(delta_w);
  if(delta_p)   free(delta_p);
  if(delta_q)   free(delta_q);
  if(delta_v)   free(delta_v);
  if(d)         free(d);
  if(e)         free(e);

  if(chol_hx)   free(chol_hx);
  if(hxinv_a)   free(hxinv_a); 
  if(P)         free(P);
}

/** Allocate internal memory dual problem related variables
 *
 */
void CL2N2_prLOQO::reallocDualMem(int mm)
{
  int nn = cur_mem_n;

  r         = (double*)realloc(r,         mm*sizeof(double));
  b_local   = (double*)realloc(b_local,   mm*sizeof(double));
  dual      = (double*)realloc(dual,      (5*mm)*sizeof(double));
  workspace = (double*)realloc(workspace, (nn*(mm+2)+2*mm)*sizeof(double));
  h_y       = (double*)realloc(h_y,       mm*mm*sizeof(double));
  c_y       = (double*)realloc(c_y,       mm*sizeof(double));
  rho       = (double*)realloc(rho,       mm*sizeof(double));
  alpha     = (double*)realloc(alpha,     mm*sizeof(double));
  beta      = (double*)realloc(beta,      mm*sizeof(double));
  gamma_w   = (double*)realloc(gamma_w,   mm*sizeof(double));
  gamma_q   = (double*)realloc(gamma_q,   mm*sizeof(double));
  hat_alpha = (double*)realloc(hat_alpha, mm*sizeof(double));
  hat_beta  = (double*)realloc(hat_beta,  mm*sizeof(double));
  delta_y   = (double*)realloc(delta_y,   mm*sizeof(double));
  delta_w   = (double*)realloc(delta_w,   mm*sizeof(double));
  delta_p   = (double*)realloc(delta_p,   mm*sizeof(double));
  delta_q   = (double*)realloc(delta_q,   mm*sizeof(double));
  delta_v   = (double*)realloc(delta_v,   mm*sizeof(double));
  e         = (double*)realloc(e,         mm*sizeof(double));

  cur_mem_m = mm;
}

/** Allocate internal memory primal problem related variables
 *
 */
void CL2N2_prLOQO::reallocPrimalMem(int nn)
{
  int mm = cur_mem_m;

  h_x       = (double*)realloc(h_x,       nn*nn*sizeof(double));
  primal    = (double*)realloc(primal,    4*nn*sizeof(double));
  workspace = (double*)realloc(workspace, (nn*(mm+2)+2*mm)*sizeof(double));
  diag_h_x  = (double*)realloc(diag_h_x,  nn*sizeof(double));
  c_x       = (double*)realloc(c_x,       nn*sizeof(double));
  h_dot_x   = (double*)realloc(h_dot_x,   nn*sizeof(double));
  nu        = (double*)realloc(nu,        nn*sizeof(double));
  tau       = (double*)realloc(tau,       nn*sizeof(double));
  sigma     = (double*)realloc(sigma,     nn*sizeof(double));
  gamma_z   = (double*)realloc(gamma_z,   nn*sizeof(double));
  gamma_s   = (double*)realloc(gamma_s,   nn*sizeof(double));
  hat_nu    = (double*)realloc(hat_nu,    nn*sizeof(double));
  hat_tau   = (double*)realloc(hat_tau,   nn*sizeof(double));
  delta_x   = (double*)realloc(delta_x,   nn*sizeof(double));
  delta_s   = (double*)realloc(delta_s,   nn*sizeof(double));
  delta_z   = (double*)realloc(delta_z,   nn*sizeof(double));
  delta_g   = (double*)realloc(delta_g,   nn*sizeof(double));
  delta_t   = (double*)realloc(delta_t,   nn*sizeof(double));
  d         = (double*)realloc(d,         nn*sizeof(double));

  chol_hx   = (double*)realloc(chol_hx,   nn*nn*sizeof(double));
  hxinv_a   = (double*)realloc(hxinv_a,   nn*sizeof(double));
  P         = (double*)realloc(P,   nn*sizeof(double));

  cur_mem_n = nn;
}



/*****************************************************************
   taken from numerical recipes and modified to accept pointers
   moreover numerical recipes code seems to be buggy (at least the
   ones on the web)

   cholesky solver and backsubstitution
   leaves upper right triangle intact (rows first order)
   ***************************************************************/
// chteo: Original choldc compute U^TU and store it in lower triangle of a

void CL2N2_prLOQO::choldc(double a[], int n, double p[])
{
#ifdef USE_LAPACK_ATLAS
        
   // make a clone of matrix a
   double* a2=new double[n*n];
   for (int i=0; i<n; i++)
   {
      for (int j=0; j<n; j++)
      {
         a2[n*i+j]=a[n*i+j];
      }
   }
   
   // cholesky decomposition
   int info=clapack_dpotrf(CblasRowMajor, CblasUpper, n, a2, n); 
   
   // store diagonal of cholesky output (upper triangle) in vector p
   for (int i=0; i<n; i++)
      p[i] = a2[(n+1)*i];
   
   // transpose the cholesky output (upper triangle) and store it in matrix a
   for (int i=0; i<n; i++)
   {
      for (int j=i+1; j<n; j++)
      {
         a[n*j+i]=a2[n*i+j];
      }
   }
   
   delete[] a2;
   
   if (info > 0)
   {
      printf("choldc(): dpotrf (Cholesky decomposition) failed with  info=%d\n",info);
//       double maxD = -1e100;
//       double minD = 1e100;
//       cout << "D:  ";
//       for(int i=0; i<n; i++)
//       {
//          maxD = max(maxD,d[i]);
//          minD = min(minD,d[i]);
//          cout << d[i] << " ";
//       }
//       cout << endl;
//       cout << "max(D) = " << maxD << endl;
//       cout << "min(D) = " << minD << endl;
      exit(EXIT_FAILURE);
   }

#else
   int i, j, k;
   double sum;
   
   for (i = 0; i < n; i++)
   {
      for (j = i; j < n; j++) 
      {
         sum=a[n*i + j];
         for (k=i-1; k>=0; k--) 
            sum -= a[n*i + k]*a[n*j + k];
         if (i == j) 
         {
            if (sum <= 0.0) 
            {        
               printf("choldc failed, matrix not positive definite!\n");
               printf("sum : %.16f\n",sum);
               exit(0);
            }
            p[i]=sqrt(sum);
         } 
         else 
            a[n*j + i] = sum/p[i];
      }
   }
#endif
}


/** solve linear system using cholesky dc output (lower triangle)
 *
 */
void CL2N2_prLOQO::cholsb(double a[], int n, double p[], double b[], double x[])
{
#ifdef USE_LAPACK_ATLAS
   // chteo: change the UPLO to lower when we change a to lower triangle, later.
   // http://atlas3.sourcearchive.com/documentation/3.6.0-20.6/clapack__dpotrs_8c-source.html
   
   for(int i=0; i<n; i++)
      a[i*n+i] = p[i];
   memcpy(x,b,sizeof(double)*n);
   int info = clapack_dpotrs(CblasRowMajor, CblasUpper, n, 1, a, n, x, n);
   
   if(info != 0)
   {
      cout << "cholsb(): dpotrs (Linear system solver with Cholesky output) failed with info = " << info << endl;
      exit(EXIT_FAILURE);
   }
   
#else
  int i, k;
  double sum;
  
  // forward pass with lower triangle
  for (i=0; i<n; i++) {
     sum=b[i];
     for (k=i-1; k>=0; k--) sum -= a[n*i + k]*x[k];
     x[i]=sum/p[i];
  }
  
  // backward pass with  upper triangle
  for (i=n-1; i>=0; i--) {
     sum=x[i];
     for (k=i+1; k<n; k++) sum -= a[n*k + i]*x[k];
     x[i]=sum/p[i];
  }
#endif
}

/*****************************************************************
  sometimes we only need the forward or backward pass of the
  backsubstitution, hence we provide these two routines separately 
  ***************************************************************/

void CL2N2_prLOQO::chol_forward(double a[], int n, double p[], double b[], double x[])
{
   int i, k;
   double sum;
   
   for (i=0; i<n; i++) {
      sum=b[i];
      for (k=i-1; k>=0; k--) sum -= a[n*i + k]*x[k];
      x[i]=sum/p[i];
   }
}

void CL2N2_prLOQO::chol_backward(double a[], int n, double p[], double b[], double x[])
{
   int i, k;
   double sum;
   
   for (i=n-1; i>=0; i--) {
      sum=b[i];
      for (k=i+1; k<n; k++) sum -= a[n*k + i]*x[k];
      x[i]=sum/p[i];
   }
}



// chteo: we specialize it to the case where there is only 1 linear constraint
// if choldc failed, try to solve reduced KKT system by lapack.dsysv

/*****************************************************************
  solves the system | -H_x A' | |x_x| = |c_x|
                    |  A   H_y| |x_y|   |c_y|

  with H_x (and H_y) positive (semidefinite) matrices
  and n, m the respective sizes of H_x and H_y

  for variables see pg. 48 of notebook or do the calculations on a
  sheet of paper again

  predictor solves the whole thing, corrector assues that H_x didn't
  change and relies on the results of the predictor. therefore do
  _not_ modify workspace

  if you want to speed tune anything in the code here's the right
  place to do so: about 95% of the time is being spent in
  here. something like an iterative refinement would be nice,
  especially when switching from double to single precision. if you
  have a fast parallel cholesky use it instead of the numrec
  implementations.

  side effects: 1. changes H_y (but this is just the unit matrix 
                   or zero anyway in our case)
                2. changes H_x (due to choldc)
  ***************************************************************/

void CL2N2_prLOQO::solve_reduced(int n, int m, double h_x[], double h_y[], 
                                 double a[], double x_x[], double x_y[],
                                 double c_x[], double c_y[],
                                 double workspace[], int step)
{
#ifdef USE_LAPACK_ATLAS
   assert(m == 1);

   int info = 0;

   // precompute and store factorization of h_x
   if(step == PRLOQO::PREDICTOR)
   {

      memcpy(chol_hx, h_x, sizeof(double)*n*n);
      info=clapack_dpotrf(CblasRowMajor, CblasUpper, n, chol_hx, n);
      //LAPACK_ERR_CHK("solve_reduced","dpotrf 1",info);
      if(info != 0)
      {
         // Jacobi preconditioning: M = PP^-1MP^-1P, P = sqrt(diag(M))
         cout << "*";
         memcpy(chol_hx, h_x, sizeof(double)*n*n);
         for(int i=0; i<n; i++)
            P[i] = sqrt(chol_hx[i*n+i]);

         for(int i=0; i<n; i++)
            for(int j=0; j<n; j++)
               chol_hx[i] = chol_hx[i*n+j] / P[i] / P[j];

         info=clapack_dpotrf(CblasRowMajor, CblasUpper, n, chol_hx, n);

         // chol_hx (upper triangle) * P
         for(int i=0; i<n; i++)
            for(int j=i; j<n; j++)
               chol_hx[i*n+j] *= P[j];
      }
      LAPACK_ERR_CHK("solve_reduced","dpotrf 1",info);

      memcpy(hxinv_a, a, sizeof(double)*n);
      info = clapack_dpotrs(CblasRowMajor, CblasUpper, n, 1, chol_hx, n, hxinv_a, n);
      LAPACK_ERR_CHK("solve_reduced","dpotrs 1",info);

      hy_aug = h_y[0] + cblas_ddot(n, hxinv_a, 1, a, 1);
   }
   
   // solve reduced KKT
   x_y[0] = (c_y[0] + cblas_ddot(n, hxinv_a, 1, c_x, 1)) / hy_aug;
   
   memcpy(x_x, c_x, sizeof(double)*n);
   info = clapack_dpotrs(CblasRowMajor, CblasUpper, n, 1, chol_hx, n, x_x, n);
   LAPACK_ERR_CHK("solve_reduced","dpotrs 2",info);
   cblas_dscal(n, -1.0, x_x, 1);
   cblas_daxpy(n, x_y[0], hxinv_a, 1, x_x, 1);

#else
   int i,j,k;
   double *p_x;
   double *p_y;
   double *t_a;
   double *t_c;
   double *t_y;
   
   p_x = workspace;		/* together n + m + n*m + n + m = n*(m+2)+2*m */
   p_y = p_x + n;
   t_a = p_y + m;
   t_c = t_a + n*m;
   t_y = t_c + n;

   assert(m == 1);

   // this is block is equivalent to the PresolveKKT() in the python loqo
   if (step == PRLOQO::PREDICTOR) 
   {
      // h_x is augmented with the diagonal matrix d before jumping into this function
      choldc(h_x, n, p_x);	/* do cholesky decomposition */
      chol_forward(h_x, n, p_x, a, t_a);
      
      // chteo: h_y is E in loqo paper
      /* compute (h_y + a h_x^-1 a') */
      for (k=0; k<n; k++) 
         h_y[0] += t_a[k] * t_a[k];
   }
   
   chol_forward(h_x, n, p_x, c_x, t_c);   /* forward pass for c */
   
   t_y[0] = c_y[0];
   for (j=0; j<n; j++)
      t_y[0] += t_a[j] * t_c[j];
   x_y[0] = c_y[0] / h_y[0]; 
   
   /* finally solve for x_x */
   for (i=0; i<n; i++)
      t_c[i] = -t_c[i] + t_a[i] * x_y[0];
   chol_backward(h_x, n, p_x, t_c, x_x);
#endif
}


/*****************************************************************
  matrix vector multiplication (symmetric matrix but only one triangle
  given). computes m*x = y
  no need to tune it as it's only of O(n^2) but cholesky is of
  O(n^3). so don't waste your time _here_ although it isn't very
  elegant. 
  ***************************************************************/

void CL2N2_prLOQO::matrix_vector(int n, double m[], double x[], double y[])
{
  int i, j;

  for (i=0; i<n; i++) {
    y[i] = m[(n+1) * i] * x[i];

    for (j=0; j<i; j++)
      y[i] += m[i + n*j] * x[j];

    for (j=i+1; j<n; j++) 
      y[i] += m[n*i + j] * x[j]; 
  }
}

/*****************************************************************
  call only this routine; this is the only one you're interested in
  for doing quadratical optimization

  the restart feature exists but it may not be of much use due to the
  fact that an initial setting, although close but not very close the
  the actual solution will result in very good starting diagnostics
  (primal and dual feasibility and small infeasibility gap) but incur
  later stalling of the optimizer afterwards as we have to enforce
  positivity of the slacks.
  ***************************************************************/
void CL2N2_prLOQO::SolveQP()
{
   int &n = dim;
   int &m = numOfConstraint;

   // check if pre allocated mem is enough
   if(n > cur_mem_n)
      reallocPrimalMem(n+gradSetSizeIncrement);
   if(m > cur_mem_m) {
      reallocDualMem(m+gradSetSizeIncrement);
   }

   double bbb = 0.0;
   double rrr = *b;  // the RHS of linear constraint

   // switch between (in)equality constraint
   if(!nonNegativeSlack)
   {
	   bbb = *b;
	   rrr = 0.0;
   }
   for(int i=0; i<m; i++) 
   {
	   b_local[i] = bbb;
	   r[i] = rrr;
   }
   

   // reset solution vector
   //memset(x, 0, sizeof(double)*n);
   
   /* the knobs to be tuned ... */
   double margin = -0.95;	/* we will go up to 95% of the
				   distance between old variables and zero */

   double bound = initBound;		/* preset value for the start. small
				   values give good initial
				   feasibility but may result in slow
				   convergence afterwards: we're too
				   close to zero */
   
   // variable name change between wrapper and original pr_loqo impl
   double *c = f;
   
   // h_x get changed after each predictor-corrector pair call
   // alternatively, copy the upper triangle back to the lower one
   memcpy(h_x, Q, n*n*sizeof(double));
   
   
   /* from the header - pointers into primal and dual */
   //double *x;  // using member variable #x#
   double *y;
   double *g;
   double *z;
   double *s;
   double *t;
   
   // for ineq constraint
   double *v;
   double *w;
   double *p;  
   double *q;


   /* auxiliary variables */
   double b_plus_1;
   double c_plus_1;
   
   double x_h_x;
   double primal_inf;
   double dual_inf;
   
   double sigfig;
   double primal_obj, dual_obj;
   double mu;
   double steplength = 0;
   int counter = 0;
   
   int status = 0;
   
   int i,j;
   
   // chteo: dont mix x with other primal variables anymore
   //x = primal;			/* n */
   //g = x + n;			/* n */
   
   
   g = primal;
   t = g + n;			/* n */
   z = t + n;			/* n */
   s = z + n;			/* n */
   
   y = dual;	            /* m */
   w = y + m;		    /* m */
   p = w + m;		    /* m */
   q = p + m;		    /* m */
   v = q + m;		    /* m */
   
   
   
   /* initial settings */
   c_plus_1 = 0;
   for(i=0; i<n; i++)
      c_plus_1 += c[i]*c[i];
   c_plus_1 = sqrt(c_plus_1) + 1;
   
   b_plus_1 = 0;
   for(i=0; i<m; i++)
      b_plus_1 += b_local[i]*b_local[i];
   b_plus_1 = sqrt(b_plus_1) + 1;
   
   
   /* get diagonal terms */
   for (i=0; i<n; i++) diag_h_x[i] = h_x[(n+1)*i]; 
   
   // chteo: no more restart, use default start settings
   for (i=0; i<m; i++)
      for (j=i; j<m; j++)
         h_y[i*m + j] = (i==j) ? 1 : 0;
   
   for (i=0; i<n; i++) {
      c_x[i] = c[i];
      h_x[(n+1)*i] += 1;
   }
   
   for (i=0; i<m; i++)
      c_y[i] = r[i];
   
   /* and solve the system [-H_x A'; A H_y] [x, y] = [c_x; c_y] */
   solve_reduced(n, m, h_x, h_y, a, x, y, c_x, c_y, workspace,
                 PRLOQO::PREDICTOR);
   
   
   /* initialize the other variables (primal) */
   for (i=0; i<n; i++) {
      g[i] = max(fabs(x[i] - l[i]), bound);
      z[i] = max(fabs(x[i]), bound);
      t[i] = max(fabs(u[i] - x[i]), bound); 
      s[i] = max(fabs(x[i]), bound); 
   }
   
   /* initialize the other variables (dual) */
   // chteo:
   // m == 1
//    for (i=0; i<m; i++) {    
//       v[i] = max(fabs(y[i]), bound);
//       w[i] = max(fabs(y[i]), bound);
//       q[i] = max(fabs(y[i]), bound); 
//       p[i] = max(fabs(r[i] - w[i]), bound);
//    }
   v[0] = max(fabs(y[0]), bound);
   w[0] = max(fabs(y[0]), bound);
   q[0] = max(fabs(y[0]), bound); 
   p[0] = max(fabs(r[0] - w[0]), bound);

   
   mu = 0;
   for(i=0; i<n; i++)
      mu += z[i] * g[i] + s[i] * t[i];
   for(i=0; i<m; i++)
      mu += v[i] * w[i] + p[i] * q[i];
   mu = mu / (2*(m+n));
   sigfig = 0;
   steplength = 1;
   
   /* the main loop */
   if (verbosity >= PRLOQO::STATUS) {
      printf("counter | pri_inf  | dual_inf | pri_obj   | dual_obj  | ");
      printf("sigfig | alpha  | nu \n");
      printf("-------------------------------------------------------");
      printf("---------------------------\n");
   }
   
   while (status == 0) {
      /* predictor */
      
      /* put back original diagonal values */
      for (i=0; i<n; i++) 
         h_x[(n+1) * i] = diag_h_x[i];

#ifdef USE_LAPACK_ATLAS
      // h_x is symmetric
      cblas_dsymv(CblasRowMajor, CblasUpper, n, 1.0, h_x, n, x, 1, 0.0, h_dot_x, 1);
      //matrix_vector(n, h_x, x, h_dot_x); /* compute h_dot_x = h_x * x */
#else
      matrix_vector(n, h_x, x, h_dot_x); /* compute h_dot_x = h_x * x */
#endif
      
      // chteo: m == 1
//       for (i=0; i<m; i++) {
//          rho[i] = b_local[i] + w[i];
//          for (j=0; j<n; j++)
//             rho[i] -= a[n*i + j] * x[j];
//       }
      rho[0] = b_local[0] + w[0];
      for (j=0; j<n; j++)
	  rho[0] -= a[j] * x[j];

      
      for (i=0; i<n; i++) {
         nu[i] = l[i] - x[i] + g[i];
         tau[i] = u[i] - x[i] - t[i];
         
         sigma[i] = c[i] - z[i] + s[i] + h_dot_x[i];

	 // m == 1
//          for (j=0; j<m; j++)
//             sigma[i] -= a[n*j + i] * y[j];
	 sigma[i] -= a[i] * y[0];
         
         gamma_z[i] = - z[i];
         gamma_s[i] = - s[i];
      }
      
      // m == 1
//       for(i=0; i<m; i++) {
//          alpha[i] = r[i] - w[i] - p[i];
//          beta[i]  = y[i] + q[i] - v[i];
         
//          gamma_w[i] = - w[i];
//          gamma_q[i] = - q[i];
//       }
      alpha[0] = r[0] - w[0] - p[0];
      beta[0]  = y[0] + q[0] - v[0];
      gamma_w[0] = - w[0];
      gamma_q[0] = - q[0];

      
      /* instrumentation */
      x_h_x = 0;
      primal_inf = 0;
      dual_inf = 0;
      
      for (i=0; i<n; i++) {
         x_h_x += h_dot_x[i] * x[i];
         primal_inf += SML::sqr(tau[i]);
         primal_inf += SML::sqr(nu[i]);
         dual_inf += SML::sqr(sigma[i]);
      }
      
      // m == 1
//       for (i=0; i<m; i++) {
//          primal_inf += SML::sqr(rho[i]);
//          primal_inf += SML::sqr(alpha[i]);
//          dual_inf += SML::sqr(beta[i]);
//       }
      primal_inf += SML::sqr(rho[0]);
      primal_inf += SML::sqr(alpha[0]);
      dual_inf += SML::sqr(beta[0]);
      
      primal_inf = sqrt(primal_inf)/b_plus_1;
      dual_inf = sqrt(dual_inf)/c_plus_1;
      
      primal_obj = 0.5 * x_h_x;
      dual_obj = -0.5 * x_h_x;
      for (i=0; i<n; i++) {
         primal_obj += c[i] * x[i];
         dual_obj += l[i] * z[i] - u[i] * s[i];
      }
      // m == 1
//       for (i=0; i<m; i++) {
//          dual_obj += b_local[i] * y[i];
//          dual_obj -= r[i] * q[i];
//       }
      dual_obj += b_local[0] * y[0];
      dual_obj -= r[0] * q[0];

      
      //sigfig = log10(fabs(primal_obj) + 1) - log10(fabs(primal_obj - dual_obj));
      //sigfig = max(sigfig, 0.0);
      sigfig = max(-log10(fabs(primal_obj-dual_obj)/(fabs(primal_obj)+1.0)), 0.0);

      /* the diagnostics - after we computed our results we will
         analyze them */
      
      if (counter >= maxIntPointIter) status = PRLOQO::ITERATION_LIMIT;
      if (sigfig  >= sigfig_max)  status = PRLOQO::OPTIMAL_SOLUTION;
      if (primal_inf >= 1e50)   status = PRLOQO::PRIMAL_INFEASIBLE;
      if (dual_inf >= 1e50)     status = PRLOQO::DUAL_INFEASIBLE;
      if ((primal_inf >= 1e50) & (dual_inf >= 1e50)) status = PRLOQO::PRIMAL_AND_DUAL_INFEASIBLE;
      if (fabs(primal_obj) >= 1e50) status = PRLOQO::PRIMAL_UNBOUNDED;
      if (fabs(dual_obj) >= 1e50) status = PRLOQO::DUAL_UNBOUNDED;
      
      /* write some nice routine to enforce the time limit if you
         _really_ want, however it's quite useless as you can compute
         the time from the maximum number of iterations as every
         iteration costs one cholesky decomposition plus a couple of 
         backsubstitutions */
      
      /* generate report */
      if ((verbosity >= PRLOQO::FLOOD) | ((verbosity == PRLOQO::STATUS) & (status != 0)))
         printf("%7i | %.2e | %.2e | % .2e | % .2e | %6.3f | %.4f | %.2e\n",
                counter, primal_inf, dual_inf, primal_obj, dual_obj,
                sigfig, steplength, mu);
      
      counter++;
      
      if (status == 0) {		
         /* 
            we may keep on going, otherwise it'll cost one loop extra plus a
            messed up main diagonal of h_x 
         */
         
         /* intermediate variables (the ones with hat) */
         for (i=0; i<n; i++) {
            hat_nu[i]    = nu[i] + g[i] * gamma_z[i] / z[i];
            hat_tau[i]   = tau[i] - t[i] * gamma_s[i] / s[i];
            /* diagonal terms */
            d[i] = z[i] / g[i] + s[i] / t[i];
         }
// 	 double maxd = -1e100;
// 	 double mind = 1e100;
// 	 for(i=0; i<n; i++)
// 	 {
// 	     maxd = max(maxd,d[i]);
// 	     mind = min(mind,d[i]);
// 	 }
// 	 cout << "maxd=" << maxd << "   mind=" << mind << endl;
 
         // m == 1
//          for(i=0; i<m; i++) {
//             hat_beta[i]  = beta[i] - v[i] * gamma_w[i] / w[i];
//             hat_alpha[i] = alpha[i] - p[i] * gamma_q[i] / q[i];
//          }
         hat_beta[0]  = beta[0] - v[0] * gamma_w[0] / w[0];
         hat_alpha[0] = alpha[0] - p[0] * gamma_q[0] / q[0];

         
         /* initialization before the cholesky solver */
         for (i=0; i<n; i++) {
            h_x[(n+1)*i] = diag_h_x[i] + d[i];
            c_x[i] = sigma[i] - z[i] * hat_nu[i] / g[i] - s[i] * hat_tau[i] / t[i];
         }           
         

         // m == 1
//          memset(h_y, 0, m*m*sizeof(double));
//          for (i=0; i<m; i++) {
//             c_y[i] = rho[i];
//             e[i] = 1 / (v[i] / w[i] + q[i] / p[i]);
//             h_y[(m+1)*i] = e[i];
//          }
         c_y[0] = rho[0];
         e[0] = 1.0 / (v[0] / w[0] + q[0] / p[0]);
         h_y[0] = e[0];

         
         /* and do it */
         solve_reduced(n, m, h_x, h_y, a, delta_x, delta_y, c_x, c_y, workspace,
                       PRLOQO::PREDICTOR);
         
         // chteo: for PRIMAL related var
         for(i=0; i<n; i++) {
            /* backsubstitution */
            delta_s[i] = s[i] * (delta_x[i] - hat_tau[i]) / t[i];
            delta_z[i] = z[i] * (hat_nu[i] - delta_x[i]) / g[i];        
            delta_g[i] = g[i] * (gamma_z[i] - delta_z[i]) / z[i];
            delta_t[i] = t[i] * (gamma_s[i] - delta_s[i]) / s[i];
            
            /* central path (corrector) */
            gamma_z[i] = mu / g[i] - z[i] - delta_z[i] * delta_g[i] / g[i];
            gamma_s[i] = mu / t[i] - s[i] - delta_s[i] * delta_t[i] / t[i];
            
            /* (some more intermediate variables) the hat variables */
            hat_nu[i] = nu[i] + g[i] * gamma_z[i] / z[i];
            hat_tau[i] = tau[i] - t[i] * gamma_s[i] / s[i];
            
            /* initialization before the cholesky */
            c_x[i] = sigma[i] - z[i] * hat_nu[i] / g[i] - s[i] * hat_tau[i] / t[i];
         }
         
         // chteo: for DUAL related var
         // m == 1
//          for(i=0; i<m; i++) {
//             /* backsubstitution */
//             delta_w[i] = - e[i] * (hat_beta[i] - q[i] * hat_alpha[i] / p[i] + delta_y[i]);
//             delta_q[i] = q[i] * (delta_w[i] - hat_alpha[i]) / p[i];
            
//             delta_v[i] = v[i] * (gamma_w[i] - delta_w[i]) / w[i];
//             delta_p[i] = p[i] * (gamma_q[i] - delta_q[i]) / q[i];
            
            
//             /* central path (corrector) */
//             gamma_w[i] = mu / v[i] - w[i] - delta_w[i] * delta_v[i] / v[i];
//             gamma_q[i] = mu / p[i] - q[i] - delta_q[i] * delta_p[i] / p[i];
            
//             /* (some more intermediate variables) the hat variables */
//             hat_alpha[i] = alpha[i] - p[i] * gamma_q[i] / q[i];
//             hat_beta[i]  = beta[i] - v[i] * gamma_w[i] / w[i];
            
//             /* initialization before the cholesky */
//             c_y[i] = rho[i] - e[i] * (hat_beta[i] - q[i] * hat_alpha[i] / p[i]);
//          }
         /* backsubstitution */
         delta_w[0] = - e[0] * (hat_beta[0] - q[0] * hat_alpha[0] / p[0] + delta_y[0]);
         delta_q[0] = q[0] * (delta_w[0] - hat_alpha[0]) / p[0];
         delta_v[0] = v[0] * (gamma_w[0] - delta_w[0]) / w[0];
         delta_p[0] = p[0] * (gamma_q[0] - delta_q[0]) / q[0];

         /* central path (corrector) */
         gamma_w[0] = mu / v[0] - w[0] - delta_w[0] * delta_v[0] / v[0];
         gamma_q[0] = mu / p[0] - q[0] - delta_q[0] * delta_p[0] / p[0];

         /* (some more intermediate variables) the hat variables */
         hat_alpha[0] = alpha[0] - p[0] * gamma_q[0] / q[0];
         hat_beta[0]  = beta[0] - v[0] * gamma_w[0] / w[0];

         /* initialization before the cholesky */
         c_y[0] = rho[0] - e[0] * (hat_beta[0] - q[0] * hat_alpha[0] / p[0]);

         
         
         /* and do it */
         solve_reduced(n, m, h_x, h_y, a, delta_x, delta_y, c_x, c_y, workspace,
                       PRLOQO::CORRECTOR);
         
         // chteo: for primal related var
         for (i=0; i<n; i++) {
            /* backsubstitution */
            delta_s[i] = s[i] * (delta_x[i] - hat_tau[i]) / t[i];
            delta_z[i] = z[i] * (hat_nu[i] - delta_x[i]) / g[i];        
            delta_g[i] = g[i] * (gamma_z[i] - delta_z[i]) / z[i];
            delta_t[i] = t[i] * (gamma_s[i] - delta_s[i]) / s[i];
         }
         
         // chteo: for dual related var
         // m == 1
//          for(i=0; i<m; i++) {
//             /* backsubstitution */
//             delta_w[i] = - e[i] * (hat_beta[i] - q[i] * hat_alpha[i] / p[i] + delta_y[i]);
//             delta_q[i] = q[i] * (delta_w[i] - hat_alpha[i]) / p[i];
//             delta_v[i] = v[i] * (gamma_w[i] - delta_w[i]) / w[i];
//             delta_p[i] = p[i] * (gamma_q[i] - delta_q[i]) / q[i];
//          }
         /* backsubstitution */
         delta_w[0] = - e[0] * (hat_beta[0] - q[0] * hat_alpha[0] / p[0] + delta_y[0]);
         delta_q[0] = q[0] * (delta_w[0] - hat_alpha[0]) / p[0];
         delta_v[0] = v[0] * (gamma_w[0] - delta_w[0]) / w[0];
         delta_p[0] = p[0] * (gamma_q[0] - delta_q[0]) / q[0];

         
         steplength = -1;
         for (i=0; i<n; i++) {
            steplength = min(steplength, delta_g[i]/g[i]);
            steplength = min(steplength, delta_t[i]/t[i]);
            steplength = min(steplength, delta_s[i]/s[i]);
            steplength = min(steplength, delta_z[i]/z[i]);
         }
         
         // m == 1
//          for (i=0; i<m; i++) {
//             steplength = min(steplength, delta_w[i]/w[i]);
//             steplength = min(steplength, delta_p[i]/p[i]);
//             steplength = min(steplength, delta_q[i]/q[i]);
//             steplength = min(steplength, delta_v[i]/v[i]);
//          }  
         steplength = min(steplength, delta_w[0]/w[0]);
         steplength = min(steplength, delta_p[0]/p[0]);
         steplength = min(steplength, delta_q[0]/q[0]);
         steplength = min(steplength, delta_v[0]/v[0]);

         steplength = margin / steplength;
         
         /* compute mu */
         mu = 0;
         for(i=0; i<n; i++)
            mu += z[i] * g[i] + s[i] * t[i];
         // m == 1
//          for(i=0; i<m; i++)
//             mu += v[i] * w[i] + p[i] * q[i];            
         mu += v[0] * w[0] + p[0] * q[0];            

         mu = mu / (2*(m+n));
         mu = mu * SML::sqr((steplength - 1) / (steplength + 10));
         
         for (i=0; i<n; i++) {
            x[i] += steplength * delta_x[i];
            g[i] += steplength * delta_g[i];
            t[i] += steplength * delta_t[i];
            z[i] += steplength * delta_z[i];
            s[i] += steplength * delta_s[i];
         }
         
         // m == 1
//          for (i=0; i<m; i++) {
//             y[i] += steplength * delta_y[i];
//             w[i] += steplength * delta_w[i];
//             p[i] += steplength * delta_p[i];
//             q[i] += steplength * delta_q[i];
//             v[i] += steplength * delta_v[i];
//          }
         y[0] += steplength * delta_y[0];
         w[0] += steplength * delta_w[0];
         p[0] += steplength * delta_p[0];
         q[0] += steplength * delta_q[0];
         v[0] += steplength * delta_v[0];
         
      }
   }
   if ((status == 1) && (verbosity >= PRLOQO::STATUS)) {
      printf("----------------------------------------------------------------------------------\n");
      printf("optimization converged\n");
   }

   
//    for(i=0; i<n; i++)
//    {
//        cout << "x=" << x[i] << "  d=" << d[i] << endl;
//    }

   
   // chteo:
   // h_x get changed after each predictor-corrector pair call
   // alternatively, copy the upper triangle back to the lower one
   memcpy(h_x, Q, n*n*sizeof(double));
}


// void CL2N2_prLOQO::PrintQP()
// {
//     ofstream Hfp("H.txt");
//     ofstream cfp("c.txt");
//     ofstream dfp("d.txt");

//     if(!Hfp.good() || !cfp.good() || !dfp.good())
//     {
//        cout << "Failed to open files to store QP problem" << endl;
//        return;
//     }
    
//     // Hessian matrix first
//     Hfp.precision(30);
//     for(int i=0; i<dim; i++)    
//     {
//        for(int j=0; j<dim; j++)
//           Hfp << Q[i*dim + j] << " ";
//        Hfp << endl;
//     }	

//     // linear part of the objective
//     cfp.precision(30);
//     for(int i=0; i<dim; i++)
//     {
//        cfp << f[i] << endl;
//     }

//     // the d vector in loqo
//     dfp.precision(30);
//     for(int i=0; i<dim; i++)
//        dfp << d[i] << endl; 
    
// }

#endif
