#include <stdio.h>
#include <string.h>
#include <math.h>

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

void gaussSeidel(double **A, double *b, double *x, int n){
  int max_iter = 10;
  double s;
  // double max_err = 1.0 + tol;
  //for(int iter = 0; iter < max_iter && max_err > tol; iter++)
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

      //norma maxima
      // double err = fabs(x[i] - old_xi);
      //norma maxima relativa
      //double err = fabs((x[i] - old_xi)/x[i])
      // if(err > max_err)
      //   max_err = err;
    }
    //erro residual
    // for(int i = 0; i < n; i++){
    //   double sum = 0.0;
    //   for(int j = 0; j < n; j++)
    //     sum += A[i][j] * x[j];
    
    //   double err = fabs(b[i] - sum);
    //   if(err > max_err)
    //     max_err = err;
    // }

    //OTIMIZAÇÕES DE sum += A[i][j] * x[j];
    // APENAS UM VETOR PARA MATRIZ sum += A[i * n + j] * x[j];
    // loop unrolling no j
    // double * restrict A > vetor A[i*n+j]
    
  }
}

void gaussSeidelTriDiag(double *dp, double *ds, double *di, double *b, double *x, int n, double tol){
  int iter_max = 25;
  double err = 0.0;
  double max_err = 1.0 + tol;
  double res = 0.0;
  
  for(int iter = 0; iter < iter_max && tol < max_err ; iter++){
    max_err = 0.0; 

    double x_old = x[0];
    x[0] = (b[0] - di[0]*x[1]) / dp[0];
    //norma maxima
    max_err = fabs(x[0] - x_old);
    //norma maxima relativa
    max_err = fabs((x[0] - x_old)/x[0]);
    for(int i = 1; i < n-1; i++){
      x_old = x[i];
      x[i] = (b[i] - di[i-1]*x[i-1] - ds[i]*x[i+1]) / dp[i];
      //norma maxima
      err = fabs(x[i] - x_old);
      //norma maxima relativa 
      err = fabs((x[i] - x_old)/x[i]);

      if(err > max_err) max_err = err;
    }
    x_old = x[n-2];
    x[n-1] = (b[n-1] - di[n-2]*x[n-2]) / dp[n-1]; 
    //norma maxima
    err = fabs(x[n-1] - x_old);
    //norma maxima relativa
    err = fabs((x[n-1] - x_old)/ x[n-1]);
    if(err > max_err) max_err = err;

    // erro residual
    for(int i = 0; i < n; i++){
      double sum = dp[i]*x[i];
      if(i > 0)  sum += di[i-1] *x[i-1];
      if(i < n-1) sum += ds[i] * x[i+1];

      res = fabs(b[i] - sum);
      if(res > max_err)
        max_err = err;
    }
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

void eliminacaoGaussTriDiag(double *di, double *ds, double *dp, double *b, int n) {
  for(int i = 0; i < n; i++){
    double m = di[i] / dp[i];
    di[i] = 0.0;
    dp[i+1] -= ds[i]*m;
    b[i+1] -= b[i]*m; 
  }
}

void eliminacaoGaussPentaDiag(double *a, double *b, double *c, double *d, double *e, double *f, int n) {
  double m = a[0]/ d[0];
  d[1] -= c[0]*m;
  c[1] -= f[0]*m;
  b[1] -= b[0]*m;

  for(int i = 1; i < n-3; i++){
    m = e[i]/a[i-1];
    a[i] -= d[i+1]*m;
    d[i+1] -= c[i-1]*m;
    c[i+1] -= f[i+1]*m;
    b[i+1] -= b[i]*m;
  }

  m = e[n-2]/a[n-3];
  a[n-2] -= d[n-2]*m;
  d[n-1] -= c[n-2]*m;
  b[n-1] -= b[n-2]*m;
}

void retrossubsTriDiag(double *di, double *ds, double *dp, double *x, double *b, int n){
  x[n-1] = b[n-1] / dp[n-1];
  for(int i = n-2; i >= 0; i--){
    x[i] = (b[i] - ds[i]*x[i+1]) / dp[i];
  }
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