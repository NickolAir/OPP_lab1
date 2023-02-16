#include <malloc.h>
#include <mpi.h>

double mult_on_vector(double *line, double *x, int N){
    int numproc, rank;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &numproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double res = 0;
    double new_x = 0;
    MPI_Scatter(line, (int) N / numproc, MPI_DOUBLE, &res, (int) N / numproc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < N; ++i) {
        res += line[i] * x[i];
    }
    MPI_Gather(&new_x, (int)N / numproc, MPI_DOUBLE, &res, (int)N / numproc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
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