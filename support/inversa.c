#include <stdio.h>
#include <string.h>
#include <math.h>

void imprimeMatriz(double **M, int n){
    for(int i = 0; i < n; i ++){
        for(int j = 0; j < n; j++){
            printf("%f ",M[i][j]);
        }
        printf("\n");
    }
}


void imprimeVetor(double *v, int n){
    for(int i = 0; i < n; i++){
        printf("%f ",v[i]);
    }
    printf("\n");
}

double **criaMatrizIdentidade(int n){
  double **I =(double **) malloc(sizeof(double *)*n);
  I[0] = malloc(sizeof(double *)*n*n);

  for(int i = 0; i < n; ++i){
     I[i] = I[0] + i*n;
  }

  for(int i=0; i< n; ++i){
    for(int j=0; j<n; ++j){
      if(i==j){
        I[i][i]=1;
      }
      else{
        I[i][j]=0;
      }
    }
  }

  return I;
}
// Correto!
double **fatoracaoLU(double **A, double **I, int n){
  for(int i =0; i<n; ++i){
    int iMax = i;
    for(int k = i+1; k<n; ++k){
      if(A[k][i] > A[iMax][i])
        iMax = k;
    }
    if(iMax != i){
      double *tmp;
      tmp = A[i];
      A[i] = A[iMax];
      A[iMax] = tmp;

      tmp = I[i];
      I[i] = I[iMax];
      I[iMax] = tmp;
    }

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

void resolveLy(double **A, double *b, double **y, int n) {
  for (int i = 0; i < n; ++i) {
    (*y)[i] = b[i];
    for (int j = 0; j < i; ++j)
      (*y)[i] -= A[i][j]*(*y)[j];
    //A[i][i] = 1
  }
}

void resolveUx(double **A, double *b, double **x, int n) {
  for (int i = n-1; i >= 0; --i) {
    (*x)[i] = b[i];

    for (int j = i+1; j < n; ++j){
      (*x)[i] -= A[i][j]* (*x)[j];
    }
    (*x)[i] /= A[i][i];  
    
  }
}


void resolveInversa(double **A, double **A_inversa, int n ){
  double **I = criaMatrizIdentidade(n);
  double **LU = fatoracaoLU(A, I, n);
  double *y, *x, *b;


  b = (double *) malloc(sizeof(double)*n);
  y = (double *) malloc(sizeof(double)*n);
  x = (double *) malloc(sizeof(double)*n);
  
  //TODO pivoteamento de I
  //pegar coluna a coluna de b
  for(int i = 0; i < n; i++){
    for(int j = 0; j< n; j++){
      // Ly = b
      b[j] = I[i][j];
    }
      resolveLy(LU, b, &y, n);
      resolveUx(LU, y, &x, n);
      for(int k = 0; k < n; k++){
        A_inversa[k][i] = x[k];
      }    
  }

}