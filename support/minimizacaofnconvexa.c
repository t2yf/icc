#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <likwid.h>
#include <math.h>

#include "utils.h"


void *alloc_vector(int size, int n){
    int total_size = size*n;
    void* ptr = malloc(total_size);

    if(!ptr)
        handle_error("Erro na alocação de memória\n");

    return ptr;
}

void **alloc_square_matrix(int size, int n){
    int total_size = size*n;
    int matriz_size = total_size*n;
    void **ptr = malloc(total_size);
    ptr[0] = malloc(matriz_size);

    for(int i = 0; i < n; ++i){
        ptr[i] = ptr[0] + i*n;
    }

    return ptr;
}

const char *optim[] = {"normal", "loop_unrolling", "gauss_seidel", "modified_newton", "one_vector_jacobian"};


// Avalia a função de Broyden F(x) no parâmetro x passado
static void calculate_Fx(real_t *F, real_t *x, int n){
    if(n < 2)
        handle_error("Tamanho n para F(x) não permitido\n");

    //f1(x)
    F[0] = x[0]*(x[0]*(-2) + 3) -(2*x[1]) + 1;

    //f2(x) - fn-1(x)
    for(int i= 1; i < n-1; i++){
        F[i] = x[i]*(x[i]*(-2) + 3) - x[i-1] - (2*x[i+1]) + 1;
    }

    //fn(x)
    F[n-1] = x[n-1]*(x[n-1]*(-2)+3) - x[n-2];
}



// Avalia a Jacobiana J no parâmetro x passado
static void calculate_Jx(real_t *J_a, real_t *J_b, real_t *J_c, real_t *x, int n){    
    for (int i = 0; i < n; ++i) J_b[i] = -4*x[i]+3;
}

static real_t * calculate_Hx(real_t *J_a, real_t *J_b, real_t *J_c, real_t *x, int n){   
    real_t *H = alloc_vector(sizeof(real_t), n); 
    for (int i = 0; i < n; ++i) H[i] = -4;
    return H;
}

// Calcula norma infinita de um vetor V tanto para ||F(xi)|| quanto para ||Delta(xi)||
// Retorna o valor máximo
static real_t calculate_maximum_norm(real_t *V, int n){
    real_t max = 0.0, abs; 
    for(int i = 0; i < n; i++){
        //TODO Otimização fabs próprio para evitar chamadas de função real_t abs_val = (V[i] < 0) ? -V[i] : V[i];
        abs = fabs(V[i]);
        if(abs > max)
            max = abs;
    }

    return max;
}

// Atualiza o vetor x com o deslocamento de delta
static void update_x(real_t *x, real_t *delta, int n){
    for(int i = 0; i < n; i++){
        x[i] += delta[i];
    }
}


// Algoritmo de eliminação de Gauss usado pelo método Eliminação de Gauss
static void gaussian_elimination(real_t *J_a, real_t *J_b, real_t *J_c, real_t *F, int n){
    for (int i = 0; i < n-1; ++i){
        #if OPTIM == ONE_VECTOR_JACOBIAN
            real_t m = -1.0/J_b[i];
            J_b[i+1] -= -2.0*m;
        #else
            real_t m = J_c[i]/J_b[i];
            J_b[i+1] -= J_a[i]*m;
        #endif

        F[i+1] -= F[i]*m;
    }
}

// Algoritmo de retro substituição usado pelo método Eliminação de Gauss
static void back_substituion(real_t *J_a, real_t *J_b, real_t *J_c, real_t *F, real_t *delta, int n){
    delta[n-1] = F[n-1]/J_b[n-1];
    for (int i = n-2; i >= 0; --i){
        #if OPTIM == ONE_VECTOR_JACOBIAN
            delta[i] = (F[i] - (-2.0 * delta[i+1]))/J_b[i];
        #else
            delta[i] = (F[i] - (J_a[i] * delta[i+1]))/J_b[i];
        #endif
    }
}

// Resolve um sistema linear Ax=b usando o método Eliminação de Gauss, adaptado para a matriz Jacobiana da tri-diagonal de Broyden
static void solve_linear_system_with_gaussian_elimination(real_t *J_a, real_t *J_b, real_t *J_c, real_t *F, real_t *delta, int n){
    gaussian_elimination(J_a, J_b, J_c, F, n);
    back_substituion(J_a, J_b, J_c, F, delta, n);
}

// Chute inicial do Delta usado no método Gauss-Seidel
static void calculate_delta_guess(real_t *delta, int n){
    for (int i = 0; i < n; ++i) delta[i] = 0.0;
}

// Resolve um sistema linear Ax=b usando o método Gauss-Seidel, adaptado para a matriz Jacobiana da tri-diagonal de Broyden
static void solve_linear_system_with_gaussian_sidel(real_t *J_a, real_t *J_b, real_t *J_c, real_t *F, real_t *delta, real_t tol, int max_iter, int n){
    // Chute inicial de Delta
    calculate_delta_guess(delta, n);

	for (int iter = 0; iter < max_iter; ++iter){
		real_t max_err = 0.0, new_delta, err_curr;
		
        // Regra geral do método
        // xi = (bi - ai1*x1 - ai2*x2 - ... - ain*xn)/aii
        // nao faz parte da soma o termo xi

        // Linha 1 não tem vizinho a esquerda
        // -4x0+3 -2 0 ... 0
        #if OPTIM == ONE_VECTOR_JACOBIAN
            new_delta = (F[0] + 2.0 * delta[1])/J_b[0];
        #else
            new_delta = (F[0] - J_a[0] * delta[1])/J_b[0];
        #endif

        err_curr = fabs(new_delta - delta[0]);
        if (err_curr > max_err) max_err = err_curr;
        delta[0] = new_delta;

        for (int i = 1; i < n-1; ++i){
            #if OPTIM == ONE_VECTOR_JACOBIAN
                new_delta = (F[i] + delta[i-1] + 2.0 * delta[i+1])/J_b[i];
            #else
                new_delta = (F[i] - J_c[i] * delta[i-1] - J_a[i] * delta[i+1])/J_b[i];
            #endif

            err_curr = fabs(new_delta - delta[i]);
            if (err_curr > max_err) max_err = err_curr;
            delta[i] = new_delta;
        }
            
        // Linha n-1 não tem vizinho a direita
        #if OPTIM == ONE_VECTOR_JACOBIAN
            new_delta = (F[n-1] + delta[n-2]) / J_b[n-1];
        #else
            new_delta = (F[n-1] - J_c[n-1] * delta[n-2]) / J_b[n-1];
        #endif

        err_curr = fabs(new_delta - delta[n-1]);
		if (err_curr > max_err) max_err = err_curr;
		delta[n-1] = new_delta;
		
		if (max_err < tol) break;
	}
}

static void opposite_Fx(real_t *Fx, real_t *b, int n){
    for(int i = 0; i < n; i++){
        b[i] = -Fx[i];
    }
}


void newton_method_min_fn_normal(int n, int max_iter, real_t epsilon, real_t *x0, FILE *output, rtime_t *newton_time, rtime_t *jacobian_time, rtime_t *sl_time){ 
    real_t *x       = alloc_vector(sizeof(real_t), n);
    real_t *delta   = alloc_vector(sizeof(real_t), n);
    real_t *Fx      = alloc_vector(sizeof(real_t), n);
    real_t *J_a     = alloc_vector(sizeof(real_t), n);
    real_t *J_b     = alloc_vector(sizeof(real_t), n);
    real_t *J_c     = alloc_vector(sizeof(real_t), n);
    real_t *b       = alloc_vector(sizeof(real_t), n); //como a eliminação de gauss altera a Fx, mandar esse vetor
    real_t *H;

    LIKWID_MARKER_START("total");
    *newton_time = timestamp();

    for(int i = 0; i < n; i++){
        x[i] = x0[i];
    }

    for(int iter = 0; iter < max_iter; iter++){
        //Derivada a primeira
        calculate_Jx(J_a, J_b, J_c, x, n);
        opposite_Fx(J_b, b, n); //b é o Fx da derivada a primeira

        //Critério de parada atendido, será retornado o x
        if (calculate_maximum_norm(J_b, n) < epsilon){
            break;
        }


        //Derivada a segunda
        H = calculate_Hx(J_a, J_b, J_c, x, n);   

        LIKWID_MARKER_START("linear_system");
        *sl_time = timestamp();
        
        solve_linear_system_with_gaussian_elimination(J_a, H, J_c, b, delta, n);
        
        *sl_time = timestamp() - *sl_time;
        LIKWID_MARKER_STOP("linear_system");

        update_x(x, delta, n);

        fprintf(output, "Iteração %d\n", iter);
        print_x(output, x, n);
        
        //Critério de parada atendido, será retornado o x novo
        if(calculate_maximum_norm(delta, n) < epsilon){
            fprintf(output, "Acabou antes\n");
            break;
        }
    }
    
    *newton_time = timestamp() - *newton_time;
    LIKWID_MARKER_STOP("total");

    free(delta);
    free(Fx);
    free(J_a);
    free(J_b);
    free(J_c);
    free(b);
}

