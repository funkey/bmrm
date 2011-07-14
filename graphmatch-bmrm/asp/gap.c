#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "gap.h"

/* based on Tiberio's matlab code
 A: size of graph G 
 I: size of graph g
*/

void print_matrix(long double **M, int A, int I){
  int a,i;
  for (a=0;a<A;a++){
    for (i=0;i<I;i++){
      printf("%d ",(int) M[a][i]);
    }
    printf("\n");
  }
  printf("\n");
}

/* since i store the matrix as a vector
 * this machinery will get the right 
 * index for me quickly */
long double get(long double *C, int A, int I, int a, int i, int b, int j){
  return C[a*A*I*I + i*A*I + b*I + j];
}

void set(long double *C, long double value, int A, int I, int a, int i, int b, int j){
  C[a*A*I*I + i*A*I + b*I + j] = value;
}

/* compute every entry of matrix Q
 *
 */
long double compute_Qai(long double *C, long double **M, int A, int I, int a, int i){
  int b,j;
  long double sum = 0.0;
  for (b=0;b<A;b++){
    for (j=0;j<I;j++){
      sum += M[b][j]*get(C,A,I,a,i,b,j);
    }
  }
  return sum;
}

/* evaluate Q
 * 
 */
void evaluate_Q(long double **Q, long double *C, long double **M, int A, int I){
  int a,i;
  for (a=0;a<A;a++){
    for (i=0;i<I;i++){
      Q[a][i] = compute_Qai(C, M, A, I, a, i);
    }
  }
  
}

/* evaluate Mt
 *
 */
void evaluate_Mt(long double **Mt, long double **Q, long double beta, int A, int I){
  int a,i;
  for (a=0;a<A;a++){
    for (i=0;i<I;i++){
      Mt[a][i] = exp(beta*Q[a][i]);
      if (beta*Q[a][i] > 50)
	Mt[a][i] = GAP_INFINITY;
    }
  }
}

/* normalize by row
 *
 */
void normalize_by_row(long double **Q, int A, int I){
  long double *v = (long double *)malloc(sizeof(long double)*A);
  int i;
  int a;
  for (a=0;a<A;a++){
    v[a] = 0;
    for (i=0;i<I;i++){
      v[a] += Q[a][i];
    }
  }
  for (a=0;a<A;a++){
    for (i=0;i<I;i++){
      if (v[a]!=0){
	Q[a][i] /= v[a];
      }
    }
  }
  free(v);
}

/* normalize by column
 *
 */
void normalize_by_column(long double **Q, int A, int I){
  long double *v = (long double *)malloc(sizeof(long double)*I);
  int a, i;  
  for (i=0;i<I;i++){
    v[i] = 0;
  }
  for (a=0;a<A;a++){
    for (i=0;i<I;i++){
      v[i] += Q[a][i];
    }
  }
  
  for (a=0;a<A;a++){
    for (i=0;i<I;i++){
      if (v[i]!=0){
	Q[a][i] /= v[i];
      }
    }
  }
  free(v);
}

/* copy Mt to Mt1
 *
 */
void copy_matrix(long double **Mt1, long double **Mt, int A, int I){
  int a, i;
  for (a=0;a<A;a++){
    for (i=0;i<I;i++){
      Mt1[a][i] = Mt[a][i];
    }
  }
}

/* compute error
 *
 */
long double compute_err(long double **Mt, long double **Mt1, int A, int I){
  int a, i;
  long double sum = 0.0;
  for (a=0;a<A;a++){
    for (i=0;i<I;i++){
      sum += fabs(Mt[a][i] - Mt1[a][i]); 
    }
  }
  return 1.0*sum/(A*I);
}

void normalize_matrix(long double *C, int A, int I, long double delta){
  int a,b;
  int i,j;
  long double sum;
  long double value;
  int total = A*A*I*I;
  long double *veryoldC = (long double *)malloc(sizeof(long double)*A*A*I*I);  
  long double *oldC     = (long double *)malloc(sizeof(long double)*A*A*I*I);
  long double error = 0;

  int iter = 0;
  do{
    /* copy the old values */
    for (i=0;i<total;i++){
      oldC[i] = C[i];
      veryoldC[i] = C[i];
    }
    for (a=0;a<A;a++){
      for (b=0;b<A;b++){
	sum = 0;
	for (i=0;i<I;i++){
	  for (j=0;j<I;j++){
	    sum += get(veryoldC,A,I,a,i,b,j);
	  }
	}
	if (sum!=0){
	  for (i=0;i<I;i++){
	    for (j=0;j<I;j++){
	      value = get(veryoldC,A,I,a,i,b,j)/sum;
	      set(oldC,value,A,I,a,i,b,j);
	    }
	  }
	}
      } 
    }
    for (i=0;i<I;i++){
      for (j=0;j<I;j++){
	sum = 0;
	for (a=0;a<A;a++){
	  for (b=0;b<A;b++){
	    sum += get(oldC,A,I,a,i,b,j);
	  }
	}
	if (sum!=0){
	  for (a=0;a<A;a++){
	    for (b=0;b<A;b++){
	      value = get(oldC,A,I,a,i,b,j)/sum;
	      set(C,value,A,I,a,i,b,j);
	    }
	  }
	}
      }
    }
    error = 0;
    for (i=0;i<total;i++){
      error += fabs(C[i] - veryoldC[i]);
    }

    iter ++;

    printf("%Lf\n", error);
    if (error < delta){
      break;
    }
  }
  while (1);
  free(oldC);
  free(veryoldC);
}
/* main routine to do graduated assignment
 *
 */
void gap(long *y, long double *C, int A, int I, long double delta){
  long double beta_0 = 0.5;
  long double beta = beta_0;
  long double beta_f = 10;
  long double beta_r = 1.075;
  int I_0 = 4;
  int I_1 = 30;
  //long double epsilon = 0.5;
  long double epsilon_0 = 0.5;
  long double epsilon_1 = 0.05;
  long double err0 = GAP_INFINITY;
  long double err1 = GAP_INFINITY;
  int i0,i1;
  long double **Mt = (long double **)malloc(sizeof(long double*)*(A+1));
  long double **Mt0 = (long double **)malloc(sizeof(long double*)*(A+1));
  long double **Mt1 = (long double **)malloc(sizeof(long double*)*(A+1));
  long double **Q = (long double **)malloc(sizeof(long double*)*A);
  int i,a;
  int maxind;
  long double maxval;


  /* normalize C first */
  //if (delta > 0) normalize_matrix(C,A,I,delta);
  
  /* allocate memory for Mt */
  for (i=0;i<A+1;i++){
    Mt[i] = (long double *)malloc(sizeof(long double)*(I+1));
    Mt0[i] = (long double *)malloc(sizeof(long double)*(I+1));
    Mt1[i] = (long double *)malloc(sizeof(long double)*(I+1));
  }
  
  /* allocate memory for Q */
  for (i=0;i<A;i++){
    Q[i] = (long double *)malloc(sizeof(long double)*I);
  }
    
  srand48(10); /* seed for random generator */
  
  /* initialize random matching matrix */
  for (a=0;a<A+1;a++){
    for (i=0;i<I+1;i++){
      Mt[a][i] = 1 + 0.25*drand48();
    }
  }
  //print_matrix(Mt,A,I);
  while (beta < beta_f){
    i0 = 0;
    while (i0 < I_0 && err0 > epsilon_0){
      i0++;
      copy_matrix(Mt0, Mt, A+1, I+1);
      evaluate_Q(Q, C, Mt, A, I);    /* evaluate Q is done without the dummy nodes */
      evaluate_Mt(Mt, Q, beta, A, I);/* evaluate Mt is done without the dummy nodes */
      i1 = 0;
      while (i1 < I_1 && err1 > epsilon_1){
	i1++;
	copy_matrix(Mt1, Mt, A+1, I+1);
	normalize_by_row(Mt, A+1, I+1);
	normalize_by_column(Mt, A+1, I+1);
	err1 = compute_err(Mt, Mt1, A+1, I+1);
      }
      err1 = GAP_INFINITY;
      err0 = compute_err(Mt, Mt0, A+1, I+1);
    }
    err0 = GAP_INFINITY;
    beta = beta*beta_r;
  }
  //printf("Mt:\n");
  //print_matrix(Mt,A,I);
  //printf("\n");
  /* try to find the matching vector */
  for (a=0;a<A;a++){
    maxval = Mt[a][0];
    maxind = 0;
    for (i=1;i<I;i++){
      if (Mt[a][i] > maxval){
	maxval = Mt[a][i];
	maxind = i;
      }
    }
    y[a] = maxind;
  }
  //printf("ybar:\n");
  //print_matrix(y, A, I);
  //printf("\n");
  /* free */
  for (i=0;i<A;i++){
    free(Q[i]);
  }
  free(Q);
  for (i=0;i<A+1;i++){
    free(Mt[i]);
    free(Mt0[i]);
    free(Mt1[i]);
  }
  free(Mt);
  free(Mt0);
  free(Mt1);
}
