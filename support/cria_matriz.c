#include <stdio.h>
#include <stdlib.h>

int main(){
    int n;
    while (scanf("%d", &n) != EOF) {
        int lin = n;
        int col = n;
        //###########FORMA DO MAZIERO
        // double **A = malloc(sizeof(double *)*lin);
        // A[0] = malloc(sizeof(double) *lin*col);
        // for(int i = 1; i < lin; i++){
        //     A[i] = A[0] + i*col;
        // }

        //##########COM APENAS UM VETOR
        double *A = (double *) malloc(n * n * sizeof(double));
        for(int i = 0; i < lin; i++){
            for(int j = 0; j < col; j++){
                scanf ("%lf", &(A[i*n+j]));
            }
        }
        // for(int i = 0; i < lin; i++){
        //     for(int j = 0; j < col; j++){
        //         scanf ("%lf", &(A[i][j]));
        //     }
        // }

        // for(int i = 0; i < lin; i++){
        //     for(int j = 0; j < col; j++){
        //         printf("%f", A[i][j]);
        //     }
        //     printf("\n");
        // }
        for(int i = 0; i < lin; i++){
            for(int j = 0; j < col; j++){
                printf("%lf", A[i*n+j]);
            }
            printf("\n");
        }
        //free(A[0]);
        free(A);
    }
    
    
    return 0;
}