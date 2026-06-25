#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"

#define DIFF 0.0

#define NRAND    ((real_t) random() / RAND_MAX)  // drand48() 
#define SRAND(a) srandom(a) // srand48(a)

real_t avaliarStyblinskiTang(real_t *x, int n) {
    real_t soma = 0.0;
    for (int i = 0; i < n; ++i) {
        soma += (x[i] * x[i] * x[i] * x[i]) - 16.0 * (x[i] * x[i]) + 5.0 * x[i];
    }
    return soma / 2.0; [cite: 7]
}

// Integral Monte Carlo da função Styblinski-Tang de 2 variáveis
real_t styblinskiTang(real_t a, real_t b, int namostras)
{
  real_t resultado;
  real_t soma = 0.0;
  
  printf("Metodo de Monte Carlo (x, y).\n");
  printf("a = (%f), b = (%f), n = (%d), variaveis = 2\n", a, b, namostras);
  
  rtime_t t_inicial = timestamp();
  
  /*
    
    AQUI IMPLEMENTE O CÁLCULO DA INTEGRAL  PELO
    MÉTODO DE MONTE CARLO
    
  */
  
  rtime_t t_final = timestamp();
  printf("Tempo decorrido: %f seg.\n", t_final - t_inicial);
  
  return resultado;
}


real_t retangulos_xy(real_t a, real_t b, int npontos) {

  real_t h;
  real_t resultado;
  real_t soma = 0.0;
  
  printf("Metodo dos Retangulos (x, y).\n");
  printf("a = (%f), b = (%f), n = (%d), h = (%lg)\n", a, b, npontos, h);
  
  rtime_t t_inicial = timestamp();

  /*  
  AQUI IMPLEMENTE O CÁLCULO DA INTEGRAL  PELO
  MÉTODO DOS RETÂNGULOS
    
  */
  int k = (int)ceil(sqrt(npontos));
  
  h = (b-a) / k;
  real_t ponto[2];

  for (int i = 0; i < k; ++i) {
    ponto[0] = a + (i + 0.5) * h; 
    
    for (int j = 0; j < k; ++j) {
      ponto[1] = a + (j + 0.5) * h; 
      
      soma += avaliarStyblinskiTang(ponto, 2);
    }
  }
  resultado = h * h * soma;
  rtime_t t_final = timestamp();
  printf("Tempo decorrido: %f seg.\n", t_final - t_inicial);
  
  return resultado;
}


int main(int argc, char **argv) {

  if (argc < 5) {
    printf("Utilização: %s inicial final n_amostras n_variaveis\n", argv[0]);
    return 1;
  }

  // INICIAR VALOR DA SEMENTE
  srandom(20252);
    
// CHAMAR FUNÇÕES DE INTEGRAÇÃO E EXIBIR RESULTADOS
  real_t res_ret = retangulos_xy(a, b, p); [cite: 9]
  printf("Resultado Retangulos (n=2): %f\n\n", res_ret);
  
  return 0;
}