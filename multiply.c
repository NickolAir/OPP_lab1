#include <malloc.h>
#include <mpi.h>

double mult_on_vector(double *line, double *x, int N){
    double res = 0;
    for (int i = 0; i < N; ++i) {
        res += line[i] * x[i];
    }
    return res;
}

double *multiply (double **A, int N, double *x){
    double *new_x = (double*) malloc(sizeof(double) * N);
    for (int i = 0; i < N; ++i) {
        new_x[i] = mult_on_vector(A[i], x, N);
    }
    return new_x;
}

void mult_on_const(double *x, int N, double t){
    for (int i = 0; i < N; ++i) {
        x[i] = x[i] * t;
    }
}