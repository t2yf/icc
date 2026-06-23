#include <stdio.h>
#include <string.h>
#include <math.h>

#include "utils.h"

void eliminacaoGauss(double **A, double *b, int n) {
  for (int i = 0; i < n; ++i) {
    int iMax = i;   	  
    for (int k = i+1; k < n; ++k) {
      if (fabs(A[k][i]) > fabs(A[iMax][i]))
	iMax = k;
    }

    if (iMax != i) {
      double *tmp, aux;
      tmp = A[i];
      A[i] = A[iMax];
      A[iMax] = tmp;

      aux = b[i];
      b[i] = b[iMax];
      b[iMax] = aux;
    }

    for (int k = i+1; k < n; ++k) {
      double m = A[k][i] / A[i][i];
      A[k][i]  = 0.0;
      for (int j = i+1; j < n; ++j)
	A[k][j] -= A[i][j]*m;
      b[k] -= b[i]*m;
    }
  }
}

void eliminacaoGaussUnrolling(double **A, double *b, int n) {
  for (int i = 0; i < n; ++i) {
    int iMax = i;   	  
    for (int k = i+1; k < n; ++k) {
      if (fabs(A[k][i]) > fabs(A[iMax][i]))
	iMax = k;
    }

    if (iMax != i) {
      double *tmp, aux;
      tmp = A[i];
      A[i] = A[iMax];
      A[iMax] = tmp;

      aux = b[i];
      b[i] = b[iMax];
      b[iMax] = aux;
    }

    for (int k = i+1; k < n - (n-i-1)%4; k+=4) {
      double aii = A[i][i];
      double m0 = A[k][i] / aii;
      double m1 = A[k+1][i] / aii;
      double m2 = A[k+2][i] / aii;
      double m3 = A[k+3][i] / aii;
      //double m4 = A[k+4][i] / aii;

      A[k][i] = 0.0;
      A[k+1][i] = 0.0;
      A[k+2][i] = 0.0;
      A[k+3][i] = 0.0;
     // A[k+4][i] = 0.0;

      for (int j = i+1; j < n; ++j)
	      A[k][j] -= A[i][j]*m0;

      for (int j = i+1; j < n; ++j)
	      A[k+1][j] -= A[i][j]*m1;
      
      for (int j = i+1; j < n; ++j)
	      A[k+2][j] -= A[i][j]*m2;
      
      for (int j = i+1; j < n; ++j)
	      A[k+3][j] -= A[i][j]*m3;
      
      // for (int j = i+1; j < n; ++j)
	    //   A[k+4][j] -= A[i][j]*m4;


      b[k] -= b[i]*m0;
      b[k+1] -= b[i]*m1;
      b[k+2] -= b[i]*m2;
      b[k+3] -= b[i]*m3;
      //b[k+4] -= b[i]*m4;
    }
    for (int k = n - (n-i-1)%4; k < n; ++k) {

      double m = A[k][i] / A[i][i];
      A[k][i]  = 0.0;
      for (int j = i+1; j < n; ++j)
	      A[k][j] -= A[i][j]*m;
      b[k] -= b[i]*m;
    }
  }

}


void eliminacaoGaussUnrollingEJam(double **A, double *b, int n) {
  for (int i = 0; i < n; ++i) {
    int iMax = i;   	  
    for (int k = i+1; k < n; ++k) {
      if (fabs(A[k][i]) > fabs(A[iMax][i]))
	iMax = k;
    }

    if (iMax != i) {
      double *tmp, aux;
      tmp = A[i];
      A[i] = A[iMax];
      A[iMax] = tmp;

      aux = b[i];
      b[i] = b[iMax];
      b[iMax] = aux;
    }

    for (int k = i+1; k < n - (n-i-1)%4; k+=4) {
      double aii =  A[i][i];
      double m0 = A[k][i] /aii;
      double m1 = A[k+1][i] / aii;
      double m2 = A[k+2][i] / aii;
      double m3 = A[k+3][i] / aii;
      //double m4 = A[k+4][i] / aii;

      A[k][i] = 0.0;
      A[k+1][i] = 0.0;
      A[k+2][i] = 0.0;
      A[k+3][i] = 0.0;
     // A[k+4][i] = 0.0;

      for (int j = i+1; j < n; ++j){
        double aij = A[i][j];
        A[k][j] -= aij*m0;
        A[k+1][j] -= aij*m1;
        A[k+2][j] -= aij*m2;
        A[k+3][j] -= aij*m3;
       // A[k+4][j] -= aij*m4;
      }
	      

      b[k] -= b[i]*m0;
      b[k+1] -= b[i]*m1;
      b[k+2] -= b[i]*m2;
      b[k+3] -= b[i]*m3;
     // b[k+4] -= b[i]*m4;
    }
    for (int k = n - (n-i-1)%4; k < n; ++k) {
      double m = A[k][i] / A[i][i];
      A[k][i]  = 0.0;
      for (int j = i+1; j < n; ++j)
	      A[k][j] -= A[i][j]*m;
      b[k] -= b[i]*m;
    }
  }
}

void gaussSeidel(double **A, double *b, double *x, int n){
  int max_iter = 10;
  double s;
  for(int iter = 0; iter < max_iter; iter++){
    printf("iter %d\n", iter);
    for(int i = 0; i < n; ++i){
      s = 0.0;
      //double old_xi = x[i];
      for(int j = 0; j < i; ++j){
        s += A[i][j] * x[j];
      }
      for(int j = i+1; j < n; ++j){
        s += A[i][j] * x[j];
      }
      x[i] = (b[i] - s) / A[i][i];

      // double err = fabs(x[i] - old_xi);
      // if(err > max_err)
      //   max_err = err;
    }
    // if(max_err < tol)
    //   break;
  }
}

void gaussJacobi(double **A, double *b, double *x, int n){
  int max_iter = 25;
  double *x_old = (double *) malloc(sizeof(double) * n);
  double s;
  for(int iter = 0; iter < max_iter; iter++){
    for(int i = 0; i < n; i++){
      x_old[i] = x[i];
    }
    for(int i = 0; i < n; ++i){
      s = 0.0;
      //double old_xi = x[i];
      for(int j = 0; j < i; ++j){
        s += A[i][j] * x_old[j];
      }
      for(int j = i+1; j < n; ++j){
        s += A[i][j] * x_old[j];
      }
      x[i] = (b[i] - s) / A[i][i];

      // double err = fabs(x[i] - old_xi);
      // if(err > max_err)
      //   max_err = err;
    }
    // if(max_err < tol)
    //   break;
  }
  free(x_old);
}


void retrossubs(double **A, double *b, double *x, int n) {
  for (int i = n-1; i >= 0; --i) {
    x[i] = b[i];
    for (int j = i+1; j < n; ++j)
      x[i] -= A[i][j]*x[j];
    x[i] /= A[i][i];
  }
}

void resolveSL(double **A, double *b, double *x, int n) {  
  eliminacaoGauss(A, b, n);
  retrossubs(A, b, x, n); 
}

void resolveSLUnrolling(double **A, double *b, double *x, int n) {  
  eliminacaoGaussUnrolling(A, b, n);
  retrossubs(A, b, x, n); 
}

void resolveSLUnrollingEJam(double **A, double *b, double *x, int n) {  
  eliminacaoGaussUnrollingEJam(A, b, n);
  retrossubs(A, b, x, n); 
}


