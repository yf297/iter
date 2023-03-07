#include <cblas.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <lapacke.h>



void chol_solve(double* A, double* b, int* n_p, int* max_iter_p, double* tol_p, double* r_norm, int* k_p){

  int n = *n_p;
  int max_iter = *max_iter_p;
  double tol = *tol_p;

  // copy of A
  double* Ac = (double*) malloc( n*n*sizeof(double));
  memcpy(Ac, A, n*n*sizeof(double));

  // 0th residual
  int k = 0;
  double* r = b;
  r_norm[k] = pow(cblas_dnrm2(n, r, 1),2);
 
  // A[0] = sqrt(A[0]);
  A[k*n + k] = pow(A[k*n + k] -  pow(cblas_dnrm2(k, A + k*n, 1), 2), 0.5);

  // first step and residual
  k = 1; 
  double x[n];
  memset(x, 0, n*sizeof(double));
  x[0] = b[0]/A[0];
  x[0] = x[0]/A[0];
  memcpy(r, b, n*sizeof(double));
  cblas_dgemv(CblasColMajor, CblasNoTrans, n, k+1, -1, Ac, n, x, 1, 1, r, 1); 
  r_norm[k] = pow(cblas_dnrm2(n, r, 1), 2); 
  
  while( (k < n) && (r_norm[k] > tol) ){
      
    // update cholesky
    int ipiv[n];
    LAPACKE_dtrtrs(LAPACK_COL_MAJOR, 'U', 'T', 'N', k, 1, A, n, A + k*n, n);
    A[k*n + k] = pow(A[k*n + k] -  pow(cblas_dnrm2(k, A + k*n, 1), 2), 0.5);
    memcpy(x, b, (k+1)*sizeof(double));
    LAPACKE_dtrtrs(LAPACK_COL_MAJOR,'U', 'T', 'N', k+1, 1, A, n, x, n);
    LAPACKE_dtrtrs(LAPACK_COL_MAJOR,'U', 'N', 'N', k+1, 1, A, n, x, n);
  
    // (k+1)st step and residual
    k = k+1;
    memcpy(r, b, n*sizeof(double));
    cblas_dgemv(CblasColMajor, CblasNoTrans, n, k+1, -1, Ac, n, x, 1, 1, r, 1); 
    r_norm[k] = pow(cblas_dnrm2(n, r, 1), 2) ; 

  }
  
  *k_p = k;
}
