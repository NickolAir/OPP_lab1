#include <mpi.h>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include "multiply.h"
#include "subtract.h"
#include <time.h>

#define E 0.00001
#define t 0.01

void print_vector (double *vect, int N){
    for (int i = 0; i < N; ++i) {
        printf("%f ", vect[i]);
    }
    printf("\n");
}

void print_matrix (double **A, int N){
    for (int i = 0; i < N; ++i) {
        print_vector(A[i], N);
    }
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

double **create_matrix (int N){
    double **A = (double **)malloc(N * sizeof(double*));
    for (int i = 0; i < N; ++i) {
        A[i] = (double*) malloc(N * sizeof(double));
    }
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j)
                A[i][j] = 2.0;
            else
                A[i][j] = 1.0;
        }
    }
    return A;
}

void delete_matrix (double **A, int N){
    for (int i = 0; i < N; ++i) {
        free(A[i]);
    }
    free(A);
}
//условие на сходимость. Если выполняется, то t со знаком +, иначе t со знаком -
int convergence (double **A, double *b, int N){
    for (int i = 0; i < N; ++i) {
        double summ = 0;
        for (int j = 0; j < N; ++j) {
            if (i != j){
                summ += A[i][j];
            }
        }
        summ += b[i];
        if (A[i][i] <= summ){
            return -1;
        }
    }
    return 1;
}

double *simple_iteration (double *Ax, double *b, int N, double *x){
    mult_on_const(Ax, N, t);
    subtract(x, Ax, N);
    return x;
}

double *multiplication (double **A, int N, double *x){
    double *new_x = (double*) malloc(sizeof(double) * N);
    for (int i = 0; i < N; ++i) {
        double tmp_sum = 0;
        for (int j = 0; j < N; ++j) {
            tmp_sum += A[i][j] * x[j];
        }
        new_x[i] = tmp_sum;
    }
    return new_x;
}

double norm (double *vector, int N){
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
    int rc;
    long double drob,drobSum = 0, Result = 0, sum;
    double startwtime = 0.0;
    double endwtime;
    int numprocs,rank;

    int N = 3;
    double **Matrix = create_matrix(N);
    double *vector = create_vector(N);
    double *x0 = create_x0(N);

    double *Ax = multiply(Matrix, N, x0);

    clock_t begin = clock();
    while (criterion(Ax, vector, N) == 0){
        x0 = simple_iteration(Ax, vector, N, x0);
        //print_vector(x0, N);
        Ax = multiply(Matrix, N, x0);
    }
    print_vector(x0, N);
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    printf("The elapsed time is %f seconds", time_spent);
    free(vector);
    free(x0);
    delete_matrix(Matrix, N);
}