#include <stdio.h>
#include <string.h>
#include <math.h>

#include "utils.h"

void eliminacaoGauss(double **A, double *b, int n) {
  for (int i = 0; i < n; ++i) { // somatorio i = 0 até n
    int iMax = i;
    for (int k = i+1; k < n; ++k)
      if (A[k][i] > A[iMax][i])
	iMax = k; // pivo não pode ser zero para não ter divisão por zero > pivoteamento parcial
    if (iMax != i) {
      double *tmp, aux;
      tmp = A[i];
      A[i] = A[iMax];
      A[iMax] = tmp;

      aux = b[i];
      b[i] = b[iMax];
      b[iMax] = aux;
    }

    for (int k = i+1; k < n; ++k) { // n-(i+1) vezes esse bloco
      double m = A[k][i] / A[i][i]; // 1 div
      A[k][i]  = 0.0;
      for (int j = i+1; j < n; ++j)
	A[k][j] -= A[i][j]*m; // n -(i+1) subtrações n-(i+1) multiplicações
      b[k] -= b[i]*m; // 1 sub e 1 multi
    } // custo de 2i² + i - 4ni + 2n² - n -1
  }
}

void retrossubs(double **A, double *b, double *x, int n) {
  for (int i = n-1; i >= 0; --i) {
    x[i] = b[i];
    for (int j = i+1; j < n; ++j)
    //n-1 subtrações n-1 multiplicações 
    //problema de cancelamento subtrativo se valores muito diferentes, overflow se multiplicação muito grande
      x[i] -= A[i][j]*x[j];
      //problema divisão por 0 se coeficiente for 0 
    x[i] /= A[i][i]; // n divisões // problemas solucionar com eq de gauss
  }
} //custo = n² - n + n = n² operações

void resolveSL(double **A, double *b, double *x, int n) {
  eliminacaoGauss(A, b, n);
  retrossubs(A, b, x, n); 
}

