#include <cblas.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void cg(double* A, double* b, int* n_p, int* max_iter_p, double* tol_p, double* r_norm, int* k_p){

  int n = *n_p;
  int max_iter = *max_iter_p;
  double tol = *tol_p;


  // 0th step and residual
  int k = 0;
  double* r = b;
  r_norm[k] = pow(cblas_dnrm2(n, r, 1),2);
 
  // get direction
  double* q  = (double *) malloc(n * sizeof(double)); 
  memcpy(q, r, n*sizeof(double));
  double* p  = (double *) malloc(n * sizeof(double)); 
  memcpy(p, q, n*sizeof(double));
  
  // get stepsize
  double* Ap = (double *) malloc(n * sizeof(double)); 
  memset(Ap, 0.0, n * sizeof(double));
  cblas_dgemv(CblasColMajor,CblasNoTrans, n, n, 1.0, A, n, p, 1, 0.0, Ap, 1);  
  double gamma = cblas_ddot(n, p, 1, Ap, 1);
  double alpha = cblas_ddot(n, p, 1, r, 1)/gamma;

  // first step and residual
  k = 1;
  cblas_daxpy(n, -1*alpha, Ap, 1, r, 1); 
  r_norm[k] = pow(cblas_dnrm2(n, r, 1), 2); 
  
  while( (k < max_iter) && (r_norm[k] > tol) ){

    // get direction
    memcpy(q, r, n*sizeof(double));
    double c = cblas_ddot(n, q, 1, Ap, 1)/gamma;
    cblas_daxpy(n, -1*c, p, 1, q, 1); 
    memcpy(p, q, n*sizeof(double));
      
    // get stepsize
    cblas_dgemv(CblasColMajor,CblasNoTrans, n, n, 1.0, A, n, p, 1, 0.0, Ap, 1); 
    gamma = cblas_ddot(n, p, 1, Ap, 1);
    alpha = cblas_ddot(n, p, 1, r, 1)/gamma;
    
    // (k+1)st step and residual
    k = k+1;
    cblas_daxpy(n, -1*alpha, Ap, 1, r, 1); 
    r_norm[k] = pow(cblas_dnrm2(n, r, 1), 2) ; 
  }
  
  *k_p = k;
  free(p);
  free(q);
  free(Ap);

}

