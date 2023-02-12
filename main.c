#include <mpi.h>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>

#define E 0.00001
#define t 0.01

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

double *simple_iteration (double **A, double *b, int N, double *x, double *new_x){
    for (int i = 0; i < N; ++i) {
        double tmp_sum = x[i];
        for (int j = 0; j < N; ++j) {
            if (i != j){
                tmp_sum -= A[i][j] * x[j] * t;
            }
        }
        tmp_sum += b[i] * t;
        new_x[i] = tmp_sum / 2.0;
    }
    return new_x;
}
//критерий окончания итераций, через эпсилон
int criterion (double **A, double *b, int N, double *x){
    for (int i = 0; i < N; ++i) {
        double tmp_sum = 0;
        for (int j = 0; j < N; ++j) {
            tmp_sum += A[i][j] * x[j];
        }
        tmp_sum -= b[i];
        tmp_sum /= b[i];
        if (tmp_sum < E)
            return 1;
    }
    return 0;
}

int main(int argc,char *argv[]){
    int rc;
    long double drob,drobSum = 0, Result = 0, sum;
    double startwtime = 0.0;
    double endwtime;
    int numprocs,rank;

    if (rc = MPI_Init(&argc, &argv))
    {
        printf("Ошибка запуска, выполнение остановлено\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs); // Получение числа процессов
    MPI_Comm_rank(MPI_COMM_WORLD,&rank); // Получение номера процесса
    if (rank == 0)
    {
        startwtime = MPI_Wtime();
    }
    printf("Hello from process #%d of %d\n",rank,numprocs);
    if (rank == 0)
    {
        printf("%Lf\n", Result);
        endwtime = MPI_Wtime();
        printf("%lf\n", (endwtime-startwtime)*1000);
    }
    MPI_Finalize(); // Завершение работы MPI
    return 0;
}