#include <cblas.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))


void Aconj(double* pmat, double* Apmat, double* gamma, double* q, int n, int k, int s, int ind){

  double c;
  for(int j=s; j < k; j++){
    c = cblas_ddot(ind, q, 1, Apmat + j*n, 1)/gamma[j];
    cblas_daxpy(n, -1*c, pmat + j*n, 1, q, 1); 
    }

}

void scg(double* A, double* b, int* n_p, int* max_iter_p, double* tol_p, double* r_norm, int* steps_p, int* prev_p){


  // dereference
  int n = *n_p;
  int max_iter = *max_iter_p;
  double tol = *tol_p;
  int steps = *steps_p;
  int prev = *prev_p;

  double* gamma  = (double*) malloc(max_iter * sizeof(double)); 
  double* Apmat  = (double*) malloc(max_iter * n * sizeof(double));
  double* pmat   = (double*) malloc(max_iter * n * sizeof(double));
  memset(Apmat, 0.0, n * sizeof(double));
  memset(pmat, 0.0, n * sizeof(double));


  int k = 0;
  int ind = 10;
  
  double* r = b;
  r_norm[k] = pow(cblas_dnrm2(n, r, 1),2);
 
  // q
  double* q  = (double *) malloc(n * sizeof(double)); 
  memset(q, 0.0, n*sizeof(double)); 
  memcpy(q, r, ind*sizeof(double));

  //memset(q + ind , 0.0, (n-ind)*sizeof(double));  
  
  // p
  double* p  = (double *) malloc(n * sizeof(double)); 
  memset(p, 0.0, n*sizeof(double)); 

  memcpy(p, q, ind*sizeof(double));
  memcpy(pmat + k*n, p, ind*sizeof(double));

  // Ap
  double* Ap = (double *) malloc(n * sizeof(double)); 
  memset(Ap, 0.0, n * sizeof(double));
  cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, 1.0, A, n, p, 1, 0.0, Ap, 1); 
  memcpy(Apmat + k*n, Ap, n*sizeof(double));

  // gamma
  gamma[k] = cblas_ddot(n, p, 1, Ap, 1);
  
  // alpha
  double alpha = cblas_ddot(n, p, 1, r, 1)/gamma[k];

  k = k + 1;
  ind  = ind + 1;
  
  // r
  cblas_daxpy(n, -1*alpha, Ap, 1, r, 1); 
  
  // ||r|| 
  r_norm[k] = pow(cblas_dnrm2(n, r, 1),2) ; 
  
  while( (k < max_iter) && (r_norm[k] > tol) ){
  
    // q
    memcpy(q, r, ind*sizeof(double));
    //memset(q + ind , 0.0, (n-ind)*sizeof(double));  

    // p
    int m = max(k-prev, 0);
    Aconj(pmat, Apmat, gamma, q, n, k, m, ind);
    memcpy(p, q, ind*sizeof(double));
    memcpy(pmat + k*n, p, ind*sizeof(double));

    // Ap
    cblas_dgemv(CblasColMajor,CblasNoTrans, n, ind, 1.0, A, n, p, 1, 0.0, Ap, 1); 
    memcpy(Apmat + k*n, Ap, n*sizeof(double));

    // gamma
    gamma[k] = cblas_ddot(n, p, 1, Ap, 1);

    // alpha
    alpha = cblas_ddot(n, p, 1, r, 1)/gamma[k];
    
    k = k + 1;
    
    // r
    cblas_daxpy(n, -1*alpha, Ap, 1, r, 1); 
 
    // ||r||
    r_norm[k] = pow(cblas_dnrm2(n, r, 1),2) ; 

   if(k%steps==0){
    ind = min(ind+10, n);
  }
  
  }

}

