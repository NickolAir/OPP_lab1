#include <mpi.h>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include "multiply.h"
#include "subtract.h"

#define E 0.00001
#define t 0.01

void print_vector (double *vect, int N){
    for (int i = 0; i < N; ++i) {
        printf("%f ", vect[i]);
    }
    printf("\n");
}

void print_matrix (double *A, int N){
    int line = 1;
    for (int i = 0; i < N * N; ++i) {
        printf("%f ", A[i]);
        if (i == line * N - 1){
            printf("\n");
            line++;
        }
    }
    printf("\n");
}

double sqrt(double n) {
    const double epsilon = 0.000001;
    double sqrt = 0, root = 0;
    while (sqrt < n) {
        root += epsilon;
        sqrt = root * root;
    }
    return sqrt;
}

double *create_vector (int N){
    double *vect = (double*) malloc(N * sizeof(double));
    for (int i = 0; i < N; ++i) {
        vect[i] = N + 1;
    }
    return vect;
}

double *create_x0 (int N){
    double *vect = (double*) malloc(N * sizeof(double));
    for (int i = 0; i < N; ++i) {
        vect[i] = 0;
    }
    return vect;
}

double *create_matrix (int N){
    int line = 0;
    double *A = (double*) malloc(N * N * sizeof(double));
    for (int i = 0; i < N * N; ++i) {
        if (i == line * N + line){
            A[i] = 2.0;
            line++;
        }else{
            A[i] = 1.0;
        }
    }
    return A;
}

//условие на сходимость. Если выполняется, то t со знаком +, иначе t со знаком -
int convergence (double *A, double *b, int N){
    int line = 0;
    for (int i = 0; i < N; ++i) {
        double summ = 0;
        for (int j = 0; j < N; ++j) {
            if (i != line * N + line){
                summ += A[i * N + j];
            }
        }
        summ += b[i];
        if (A[line * N + line] <= summ){
            return -1;
        }
    }
    return 1;
}

double *simple_iteration (double *Ax, int N, double *x){
    mult_on_const(Ax, N, t);
    subtract(x, Ax, N);
    return x;
}

double norm (const double *vector, int N){
    double res = 0;
    for (int i = 0; i < N; ++i) {
        res += vector[i] * vector[i];
    }
    return sqrt(res);
}

//критерий окончания итераций, через эпсилон
int criterion (double *Ax, double *b, int N){
    subtract(Ax, b, N);
    double result = norm(Ax, N) / norm(b, N);
    if (result < E){
        return 1;
    } else {
        return 0;
    }
}

int main(int argc,char *argv[]){
    double start_time = 0.0;
    double end_time;
    int numprocs, rank;

    int N = 3;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int range = N / numprocs;

//    double *Ax = multiply(Matrix, N, x0);
//
//    while (criterion(Ax, vector, N) == 0){
//        x0 = simple_iteration(Ax, N, x0);
//        Ax = multiply(Matrix, N, x0);
//    }

    double *Matrix = create_matrix(N);
    double *vector_b = create_vector(N);
    double *x0 = create_x0(N);

    double *Ax = mult_on_vector(Matrix, x0, N, rank, numprocs);
    criterion(Ax, vector_b, N);
    x0 = simple_iteration(Ax, N, x0);
    Ax = mult_on_vector(Matrix, x0, N, rank, numprocs);

    criterion(Ax, vector_b, N);
    x0 = simple_iteration(Ax, N, x0);
    Ax = mult_on_vector(Matrix, x0, N, rank, numprocs);

    if (rank == 0){
        print_vector(x0, N);
        free(Matrix);
        free(vector_b);
        free(x0);
    }
    MPI_Finalize();
}