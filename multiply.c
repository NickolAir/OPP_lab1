#include <malloc.h>
#include <mpi.h>

#define ROOT 0

void mult_on_vector(double *line, double *x, int N, int rank, int numprocs, double *res){
    double sum;
    int range = N / numprocs;
    double *tmp_res = (double*) malloc(range * sizeof(double));
    double *tmp_line = (double*) malloc(N * range * sizeof(double));

    MPI_Scatter(line, range * N, MPI_DOUBLE, tmp_line, range * N, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

    for (int i = 0; i < range; ++i) {
        sum = 0;
        for (int j = 0; j < N; ++j) {
            sum += tmp_line[i * N + j] * x[j];
        }
        tmp_res[i] = sum;
    }

    MPI_Gather(tmp_res, range, MPI_DOUBLE, res, range, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

    for (int i = range * numprocs; i < N; ++i) {
        if (rank == 0 && i % numprocs != 0){
            MPI_Recv(&res[i], 1, MPI_DOUBLE, MPI_ANY_SOURCE, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank == i % numprocs){
            sum = 0;
            for (int j = 0; j < N; ++j) {
                sum += line[i * N + j] * x[j];
            }
            if (rank != 0){
                MPI_Send(&sum, 1, MPI_DOUBLE, ROOT, 123, MPI_COMM_WORLD);
            } else {
                res[i] = sum;
            }
        }
    }
    MPI_Bcast(res, N, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    free(tmp_res);
    free(tmp_line);
}

void mult_on_const(double *x, int N, double t){
    for (int i = 0; i < N; ++i) {
        x[i] = x[i] * t;
    }
}