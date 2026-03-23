#include <stdio.h>
#include <stdlib.h>
#include "resolveinversa.h"


int main(){
    int n;
    printf("Tamanho da matriz:\n");
    scanf("%d", &n);

    double **A = (double **) malloc(sizeof(double*)*n); 
    double **A_inversa = (double **) malloc(sizeof(double*)*n); 
    for(int i = 0; i < n; i++){
        A[i] = (double *) malloc(sizeof(double)*n);
        A_inversa[i] = (double *) malloc(sizeof(double)*n);
    }

    // ler matriz A
    printf("Entre com a matriz\n");
    for(int i = 0; i < n; i ++){
        for(int j = 0; j < n; j++){
            scanf("%lf", &A[i][j]);
        }
    }

    printf("Matriz A\n");
    imprimeMatriz(A, n);
    resolveInversa(A, A_inversa, n);
    printf("Matriz A⁻¹\n");
    imprimeMatriz(A_inversa, n);
    
    return 0;
}

