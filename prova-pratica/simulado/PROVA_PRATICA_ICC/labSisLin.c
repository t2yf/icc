#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "linear.h"
#include <likwid.h>

int main ()
{
  double **A, *x, *b; // sistema linear Ax = b de tamanho n
  int n, i = 1;
  int count = 0;
  LIKWID_MARKER_INIT;
  while (scanf("%d", &n) != EOF) {
  
    // aloca e lê sistema linear
    double erro;
    scanf("%lf", &erro); // nao usado
    
    A = (double **) malloc(sizeof(double *) * n);
    for (int i = 0; i < n; ++i)
      A[i] = (double *) malloc(sizeof(double) * n);
    b = (double *) malloc(sizeof(double) * n);
    x = (double *) malloc(sizeof(double) * n);
        
    for(int i = 0; i < n; ++i)
      for(int j = 0; j < n; ++j)
	scanf ("%lg", &(A[i][j]));
    for(int i=0; i < n; ++i)
      scanf ("%lg", &(b[i]));
    
    #if V == 0
      // resolve SL
      LIKWID_MARKER_START("normal");
      rtime_t ini = timestamp();
      resolveSL(A, b, x, n);
      ini = timestamp() - ini;
      LIKWID_MARKER_STOP("normal");
      FILE *arq;
      arq = fopen("time.csv", "a+");
      fprintf(arq, "%d, %s, %lf\n", count, "normal", ini);
      fclose(arq);
    
    #elif V == 1
      // resolve SL
      LIKWID_MARKER_START("unrolling");
      rtime_t ini = timestamp();
      resolveSLUnrolling(A, b, x, n);
      ini = timestamp() - ini;
      LIKWID_MARKER_STOP("unrolling");
      FILE *arq;
      arq = fopen("time.csv", "a+");
      fprintf(arq, "%s, %lf\n", "unrolling", ini);
      fclose(arq);
      printf("%d", count);

    #elif V == 2
      // resolve SL
      LIKWID_MARKER_START("unrollingEJam");
      rtime_t ini = timestamp();
      resolveSLUnrollingEJam(A, b, x, n);
      ini = timestamp() - ini;
      LIKWID_MARKER_STOP("unrollingEJam");
      FILE *arq;
      arq = fopen("time.csv", "a+");
      fprintf(arq, "%s, %lf\n", "unrollingEJam", ini);
      fclose(arq);
    
    #elif V == 3
      // resolve SL
      LIKWID_MARKER_START("gaussSeidel");
      rtime_t ini = timestamp();
      printf("gauss seidel\n");
      gaussSeidel(A, b, x, n);
      ini = timestamp() - ini;
      LIKWID_MARKER_STOP("gaussSeidel");
      FILE *arq;
      arq = fopen("time.csv", "a+");
      fprintf(arq, "%s, %lf\n", "gaussSeidel", ini);
      fclose(arq);
    
    #elif V == 4
      // resolve SL
      LIKWID_MARKER_START("gaussJacobi");
      rtime_t ini = timestamp();
      gaussJacobi(A, b, x, n);
      ini = timestamp() - ini;
      LIKWID_MARKER_STOP("gaussJacobi");
      FILE *arq;
      arq = fopen("time.csv", "a+");
      fprintf(arq, "%s, %lf\n", "gaussJacobi", ini);
      fclose(arq);
    
    #endif
    
    // imprime solucao
    printf("\n***** Sistema %d --> n = %d\n", i++, n);
    printf("  --> X: ");
    for (int j = 0; j < n; ++j)
      printf("%.7g ", x[j]);
    printf("\n");
  
    // libera SL
    free(x); 
    free(b);
    for (int j = 0; j < n; ++j)
      free(A[j]);
    free(A);
    count++;
  }
  LIKWID_MARKER_CLOSE;
  
  return 0;  
}

