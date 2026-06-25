#include <stdio.h>
#include <string.h>
#include <math.h>

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


void gaussJacobiUnrollEJam(double **A, double *b, double *x, int n){
  int max_iter = 10;
  double s, s1, s2, s3;
  double *x_old = (double *) malloc(sizeof(double) * n);
  for(int iter = 0; iter < max_iter; iter++){

    for(int i = 0; i < n-n%4; i+=4){
      x_old[i] = x[i];
      x_old[i+1] = x[i+1];
      x_old[i+2] = x[i+2];
      x_old[i+3] = x[i+3];
    }for(int i = n-n%4; i < n; i++){
      x_old[i] = x[i];
    }

    for(int i = 0; i < n-n%4; i+=4){
      s = 0.0;
      s1 = 0.0;
      s2 = 0.0;
      s3 = 0.0;

      for(int j = 0; j < i; ++j){
        double xj = x_old[j];
        s += A[i][j] * xj;
        s1 += A[i+1][j] * xj;
        s2 += A[i+2][j] * xj;
        s3 += A[i+3][j] * xj;
      }
      for(int j = i+4; j < n; ++j){
        double xj = x_old[j];
        s += A[i][j] * xj;
        s1 += A[i+1][j] * xj;
        s2 += A[i+2][j] * xj;
        s3 += A[i+3][j] * xj;
      }

      //ajuste
      double xj0 = x_old[i];
      // s0 pula xj0 (diagonal)
      s1 += A[i+1][i] * xj0;
      s2 += A[i+2][i] * xj0;
      s3 += A[i+3][i] * xj0;

      double xj1 = x_old[i+1];
      s += A[i][i+1] * xj1;
      // s1 pula xj1 (diagonal)
      s2 += A[i+2][i+1] * xj1;
      s3 += A[i+3][i+1] * xj1;

      double xj2 = x_old[i+2];
      s += A[i][i+2] * xj2;
      s1 += A[i+1][i+2] * xj2;
      // s2 pula xj2 (diagonal)
      s3 += A[i+3][i+2] * xj2;

      double xj3 = x_old[i+3];
      s += A[i][i+3] * xj3;
      s1 += A[i+1][i+3] * xj3;
      s2 += A[i+2][i+3] * xj3;
      // s3 pula xj3 (diagonal)

      x[i] = (b[i] - s) / A[i][i];
      x[i+1] = (b[i+1] - s1) / A[i+1][i+1];
      x[i+2] = (b[i+2] - s2) / A[i+2][i+2];
      x[i+3] = (b[i+3] - s3) / A[i+3][i+3];
    }
    for(int i = n-n%4; i < n; i++){
        s = 0.0;
        for(int j = 0; j < i; ++j){
          s += A[i][j] * x_old[j];
        }
        for(int j = i+1; j < n; ++j){
          s += A[i][j] * x_old[j];
        }
        x[i] = (b[i] - s) / A[i][i];
    }
  }
  free(x_old);
}


void resolveSLUnrolling(double **A, double *b, double *x, int n) {  
  eliminacaoGaussUnrolling(A, b, n);
  retrossubs(A, b, x, n); 
}

void resolveSLUnrollingEJam(double **A, double *b, double *x, int n) {  
  eliminacaoGaussUnrollingEJam(A, b, n);
  retrossubs(A, b, x, n); 
}
