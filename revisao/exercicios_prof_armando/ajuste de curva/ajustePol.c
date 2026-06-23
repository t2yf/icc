#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <fenv.h>
#include <math.h>
#include <stdint.h>

#include <likwid.h>

#include "utils.h"

/////////////////////////////////////////////////////////////////////////////////////
//   AJUSTE DE CURVAS
/////////////////////////////////////////////////////////////////////////////////////

void montaSL(double **A, double *b, int n, long long int p, double *xy) {
    
    for(int i = 0; i < n; i++){
        b[i] = 0.0;
        for(int j=0; j<n; j++){
            A[i][j] = 0.0;
        }
    }

    A[0][0] = (double)p;

    for(int i = 0; i < p; ++i){
        double x = xy[2*i];
        double y = xy[2*i+1];

        b[0] += y;
        double pow_x = 1.0;

        for(int j = 1; j < n; ++j){
            pow_x *= x;
            A[0][j] += pow_x;
            b[j] += pow_x*y;
        }

        for(int k = 1; k<n; ++k){
            pow_x *= x;
            A[k][n-1] += pow_x;
        }

    }

    for(int i = 1; i < n; i++){
        for(int j = 0; j < n-1; j++){
            A[i][j] = A[i-1][j+1];
        }
    }
//   for (int i = 0; i < n; ++i)
//     for (int j = 0; j < n; ++j) {
//       A[i][j] = 0.0;
//       for (long long int k = 0; k < p; ++k) {
// 	A[i][j] += pow(x[k], i+j);
//       }
//     }

//   for (int i = 0; i < n; ++i) {
//     b[i] = 0.0;
//     for (long long int k = 0; k < p; ++k)
//       b[i] += pow(x[k],i) * y[k];
//   }
}

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

    for (int k = i+1; k < n-(n-i-1)%4; k+=4) {
        double aii = A[i][i];
        double m0 = A[k][i] / aii;
        double m1 = A[k+1][i] / aii;
        double m2 = A[k+2][i] / aii;
        double m3 = A[k+3][i] / aii;

        A[k][i] = 0.0;
        A[k+1][i] = 0.0;
        A[k+2][i] = 0.0;
        A[k+3][i] = 0.0;

        for(int j = i+1; j < n; j++){
            double aij = A[i][j];
            A[k][i] -= aij*m0;
            A[k+1][i] -= aij*m1;
            A[k+2][i] -= aij*m2;
            A[k+3][i] -= aij*m3;
        }

        b[k] -= b[i]*m0;
        b[k+1] -= b[i]*m0;
        b[k+2] -= b[i]*m0;
        b[k+3] -= b[i]*m0;
    } for(int k = n-(n-i-1)%4; k<n; k++){
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

double Pol(double x, int G, double *alpha) {
  double Px = alpha[0];
  for (int i = 1; i <= G; ++i)
    Px += alpha[i]*pow(x,i);
  
  return Px;
}

int main() {

  int G, g; // G -> grau do polinomio
  long long int P, p; // P -> no. de pontos
  string_t marker;

  scanf("%d %lld", &G, &P);
  p = P;   // quantidade de pontos
  g = G+1; // tamanho do SL (G + 1)

  //Otimização
  double *xy = (double *) malloc(sizeof(double)*2*p);

  //double *x = (double *) malloc(sizeof(double)*p);
  //double *y = (double *) malloc(sizeof(double)*p);

  // ler numeros
  //Otimização
  for (long long int i = 0; i < p; ++i)
    scanf("%lf %lf", xy+(2*i), xy+(2*i+1));

  double **A = (double **) malloc(sizeof(double *)*g);
  for (int i = 0; i < g; ++i)
    A[i] = (double *) malloc(sizeof(double)*g);
  
  double *b = (double *) malloc(sizeof(double)*g);
  double *alpha = (double *) malloc(sizeof(double)*g); // coeficientes ajuste

  LIKWID_MARKER_INIT;
  
  // (A) Gera SL
  marker = markerName("SL",p);
  LIKWID_MARKER_START (marker);
  double tSL = timestamp();
  montaSL(A, b, g, p, xy);
  tSL = timestamp() - tSL;
  LIKWID_MARKER_STOP(marker);
  free(marker);

  // (B) Resolve SL
  marker = markerName("EG",p);
  LIKWID_MARKER_START(marker);
  double tEG = timestamp();
  eliminacaoGauss(A, b, g); 
  retrossubs(A, b, alpha, g); 
  tEG = timestamp() - tEG;
  LIKWID_MARKER_STOP(marker);
  free(marker);

  LIKWID_MARKER_CLOSE;

  // Imprime coeficientes
  for (int i = 0; i < g; ++i)
    printf("%1.15e ", alpha[i]);
  puts("");

  // Imprime resíduos
  for (long long int i = 0; i < p; ++i)
    printf("%1.15e ", fabs(xy[2*i+1] - Pol(xy[2*i],G,alpha)) );
  puts("");

  // Imprime os tempos
  printf("%lld %1.10e %1.10e\n", P, tSL, tEG);

  return 0;
}