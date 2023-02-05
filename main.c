#include <mpi.h>
#include <stdio.h>
#include <malloc.h>

double **create_matrix (int N){
    double **A = (double **)malloc(N * sizeof(double*));
    for (int i = 0; i < N; ++i) {
        A[i] = (double*) malloc(N * sizeof(double));
    }
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if(i == j)
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

int main(int argc, char** argv) {


    /*// Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    printf("Hello world from processor %s, rank %d out of %d processors\n",
           processor_name, world_rank, world_size);

    // Finalize the MPI environment.
    MPI_Finalize();*/
}