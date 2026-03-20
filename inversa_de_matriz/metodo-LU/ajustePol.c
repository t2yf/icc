#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <fenv.h>
#include <math.h>
#include <stdint.h>

#include "linear.h"
#include "utils.h"

/////////////////////////////////////////////////////////////////////////////////////
//   AJUSTE DE CURVAS
/////////////////////////////////////////////////////////////////////////////////////

void montaSL(double **A, double *b, int n, long long int p, double *x, double *y) {
    //inicializar
    for (int i = 0; i < n; ++i){
      b[i] = 0.0;
      for (int j = 0; j < n; ++j) 
      A[i][j] = 0.0;
    }
    //o primeiro somatório eh o numero de pontos
    A[0][0] = (double)p;

    //passar por todos os pontos da amostra
    for(int i = 0; i < p; ++i){
      double pow_x;
      //somar todos os y[i] (multiplicado por 1)
      b[0] += y[i];

      pow_x = 1.0;
      //preencher primeira linha e a coluna de b
      for(int j = 1; j < n; ++j){ //passar pelo grau
        pow_x *= x[i];
        A[0][j] += pow_x;
        b[j] += pow_x*y[i];
      }

      //preencher ultima coluna
      for(int k = 1; k < n; ++k){
        pow_x *= x[i];
        A[k][n-1] += pow_x;
      }
    }
    //preencher restante da matriz
    for(int i = 1; i < n; ++i){
      for(int j = 0; j < n-1; ++j){
        A[i][j] = A[i-1][j+1];
      }
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

  scanf("%d %lld", &G, &P);
  p = P;   // quantidade de pontos
  g = G+1; // tamanho do SL (G + 1)

  double *x = (double *) malloc(sizeof(double)*p);
  double *y = (double *) malloc(sizeof(double)*p);

  // ler numeros
  for (long long int i = 0; i < p; ++i)
    scanf("%lf %lf", x+i, y+i);

  double **A = (double **) malloc(sizeof(double *)*g);
  for (int i = 0; i < g; ++i)
    A[i] = (double *) malloc(sizeof(double)*g);
  
  double *b = (double *) malloc(sizeof(double)*g);
  double *alpha = (double *) malloc(sizeof(double)*g); // coeficientes ajuste
  
  // (A) Gera SL
  double tSL = timestamp();
  montaSL(A, b, g, p, x, y);
  tSL = timestamp() - tSL;

  // (B) Resolve SL
  double tEG = timestamp();
  resolveSL(A, b, alpha, g);
  tEG = timestamp() - tEG;

  // Imprime coeficientes
  printf("Coeficientes\n");
  for (int i = 0; i < g; ++i)
    printf("%1.15e ", alpha[i]);
  puts("");

  // Imprime resíduos
  printf("Resíduos\n");
  for (long long int i = 0; i < p; ++i)
    printf("%1.15e ", fabs(y[i] - Pol(x[i],G,alpha)) );
  puts("");

  // Imprime os tempos
  printf("Número de pontos: %lld Tempo de Montagem: %1.10e Tempo de resolução: %1.10e\n", P, tSL, tEG);

  return 0;
}
