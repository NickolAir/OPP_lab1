#include <malloc.h>
#include <mpi.h>

double *mult_on_vector(double *line, double *x, int N, int rank, int numprocs){
    double sum;
    double *tmp_line = (double*) malloc(sizeof(double) * N);
    double *new_x = (double*) malloc(sizeof(double) * N );
    int range = N / numprocs;
    double *tmp_res = (double*) malloc(range * sizeof(double));

    MPI_Scatter(line, range * N, MPI_DOUBLE, tmp_line, range * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i = 0; i < range; ++i) {
        sum = 0;
        for (int j = 0; j < N; ++j) {
            sum += tmp_line[i * N + j] * x[j];
        }
        tmp_res[i] = sum;
    }

    MPI_Gather(tmp_res, range, MPI_DOUBLE, new_x, range, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i = range * numprocs; i < N; ++i) {
        if (rank == 0 && i % numprocs != 0){
            MPI_Recv(&new_x[i], 1, MPI_DOUBLE, MPI_ANY_SOURCE, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank == i % numprocs){
            sum = 0;
            for (int j = 0; j < N; ++j) {
                sum += line[i * N + j] * x[j];
            }
            if (rank != 0){
                MPI_Send(&sum, 1, MPI_DOUBLE, 0, 123, MPI_COMM_WORLD);
            } else {
                new_x[i] = sum;
            }
        }
    }
    MPI_Bcast(new_x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    free(tmp_res);
    free(tmp_line);
    return new_x;
}

/*double *multiply (double **A, int N, double *x, int rank, int numprocs){
    double *new_x = (double*) malloc(sizeof(double) * N);
    for (int i = 0; i < N; ++i) {
        new_x[i] = mult_on_vector(A[i], x, N, rank, numprocs);
    }
    return new_x;
}*/

void mult_on_const(double *x, int N, double t){
    for (int i = 0; i < N; ++i) {
        x[i] = x[i] * t;
    }
}