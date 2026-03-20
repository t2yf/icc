#include <stdio.h>
#include <string.h>
#include <math.h>

#include "utils.h"

void eliminacaoGauss(double **A, double *b, int n) {
  for (int i = 0; i < n; ++i) {
    int iMax = i;
    for (int k = i+1; k < n; ++k)
      if (A[k][i] > A[iMax][i])
	iMax = k;
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

double **criaMatrizIdentidade(int n){
  double **I =(double **) malloc(sizeof(double *)*n);
  I[0] = malloc(sizeof(double *)*n*n);

  for(int i = 0; i < n; ++i){
     I[i] = I[0] + i*n;
  }

  for(int i=0; i< n; ++i){
    for(int j=0; j<n; ++j){
      if(i=j){
        I[i][i]=1;
      }
      else{
        I[i][j]=0;
      }
    }
  }

  return I;
}


double **fatoracaoLU(double **A, int n){
  double **I = criaMatrizIdentidade(n);

  for(int i =0; i<n; ++i){
    int iMax = i;
    //TODO fazer pivoteamento parcial
    // for(int k = i+1; k<n; ++k){
    //   if(A[k][i] > A[iMax][i])
    //     iMax = k;
    // }
    // if(iMax != i){
    //   //TODO melhoria só utilizar um vetor que aponta para as trocas 
    //   double *tmp, aux;
    //   tmp = A[i];
    //   A[i] = A[iMax];
    //   A[iMax] = tmp;

    //   tmp = I[i];
    //   I[i] = I[iMax];
    //   I[iMax] = tmp;
    // }

    for(int k = i+1; k < n; ++k){
      double m = A[k][i] / A[i][i];
      //atribuir multiplicador -> parte U
      A[k][i] = m;
      //escalonar -> parte L
      for (int j = i+1; j < n; ++j){
	      A[k][j] -= A[i][j]*m;
      }
    }
  }
  return A;
}

void resolveInversa(double **A, double **I, int n ){
  double **LU = fatoracaoLU(A, n);
  double *y;
  //criar vetor b

   
  //criar vetor y

  //criar vetor x

  //Ux
  //percorrer U e percorrer x
  for()

}

