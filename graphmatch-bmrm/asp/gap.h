#include <stdlib.h>
#include <math.h>
#define GAP_INFINITY 999999999999999.99
void print_matrix(long double **M, int A, int I);
long double get(long double *C, int A, int I, int a, int i, int b, int j);
long double compute_Qai(long double *C, long double **M, int A, int I, int a, int i);
void evaluate_Q(long double **Q, long double *C, long double **M, int A, int I);
void evaluate_Mt(long double **Q, long double **Mt, long double beta, int A, int I);
void normalize_matrix(long double *C, int A, int I, long double delta);
void normalize_by_row(long double **Q, int A, int I);
void normalize_by_column(long double **Q, int A, int I);
void copy_matrix(long double **Mt1, long double **Mt, int A, int I);
long double compute_err(long double **Mt, long double **Mt1, int A, int I);
void gap(long *y, long double *C, int A, int I, long double delta);
