#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define E 0.00001
#define t 0.00001
#define ROOT 0
#define DEFAULT 3

void print_vector (double *vect, int N){
    for (int i = 0; i < N; ++i) {
        printf("%f ", vect[i]);
    }
    printf("\n");
}

double *create_vectorB (int N){
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

double norm (const double *vector, int N){
    double res = 0;
    for (int i = 0; i < N; ++i) {
        res += vector[i] * vector[i];
    }
    return sqrt(res);
}

//критерий окончания итераций через эпсилон
int criterion (double *Ax, double normB, int N){
    double result = norm(Ax, N) / normB;
    if (result < E){
        return 1;
    } else {
        return 0;
    }
}

void FreeProcess(double* Matrix, double* vectorB, double* x0, double* res, double* proc_res, double* proc_line, double* proc_x, int rank) {
    if (rank == 0)
        free(Matrix);
    free(vectorB);
    free(x0);
    free(res);
    free(proc_line);
    free(proc_res);
    free(proc_x);
}

void DataDistribution(double* Matrix, double* x0, double* proc_line, int size, int numRow, int numprocs, int rank){
    int *SendNum;
    int *SendIndx;
    int RestLines = size;
    MPI_Bcast(x0, size, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

    SendIndx = (int*) malloc(numprocs * sizeof(int));
    SendNum = (int*) malloc(numprocs * sizeof(int));
    numRow = size / numprocs;
    SendNum[0] = numRow * size;
    SendIndx[0] = 0;
    for (int i = 1; i < numprocs; ++i) {
        RestLines -= numRow;
        numRow = RestLines / (numprocs - i);
        SendNum[i] = numRow * size;
        SendIndx[i] = SendIndx[i-1] + SendNum[i-1];
    }
    MPI_Scatterv(Matrix, SendNum, SendIndx, MPI_DOUBLE, proc_line, SendNum[rank],
                MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    free(SendNum);
    free(SendIndx);
}

void ParallelCalculate(double* proc_line, double* x0, double* proc_res, double* vectorB, int size, int numRow, int numprocs, int rank){
    for (int i = 0; i < numRow; i++) {
        proc_res[i] = 0 - vectorB[i];
        for (int j = 0; j < size; j++) {
            proc_res[i] += proc_line[i * size + j] * x0[j];
        }
    }
}

void ResultReply(double* proc_res, double* res, int size, int numprocs, int rank){
    int *RecvNum;
    int *RecvInd;
    int RestLines = size;

    RecvNum = (int*) malloc(numprocs * sizeof(int));
    RecvInd = (int*) malloc(numprocs * sizeof(int));

    RecvInd[0] = 0;
    RecvNum[0] = size / numprocs;
    for (int i = 1; i < numprocs; ++i) {
        RestLines -= RecvNum[i-1];
        RecvNum[i] = RestLines / (numprocs - i);
        RecvInd[i] = RecvInd[i-1] + RecvNum[i-1];
    }
    MPI_Allgatherv(proc_res, RecvNum[rank], MPI_DOUBLE, res, RecvNum, RecvInd,
                  MPI_DOUBLE, MPI_COMM_WORLD);
    free(RecvNum);
    free(RecvInd);
}

void DataDistributionSub(double* res, double* x0, double* proc_x, double* proc_res, int size, int numRow, int numprocs, int rank){
    int *SendNum;
    int *SendIndx;
    int RestLines = size;
    MPI_Bcast(x0, size, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

    SendIndx = (int*) malloc(numprocs * sizeof(int));
    SendNum = (int*) malloc(numprocs * sizeof(int));
    numRow = size / numprocs;
    SendNum[0] = numRow;
    SendIndx[0] = 0;
    for (int i = 1; i < numprocs; ++i) {
        RestLines -= numRow;
        numRow = RestLines / (numprocs - i);
        SendNum[i] = numRow;
        SendIndx[i] = SendIndx[i-1] + SendNum[i-1];
    }
    MPI_Scatterv(res, SendNum, SendIndx, MPI_DOUBLE, proc_res, SendNum[rank],
                 MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Scatterv(x0, SendNum, SendIndx, MPI_DOUBLE, proc_x, SendNum[rank],
                 MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    free(SendNum);
    free(SendIndx);
}

void ParallelSubtract(double* proc_x, double* proc_res, int size, int numRow, int numprocs, int rank){
    for (int i = 0; i < numRow; i++) {
        proc_x[i] -= proc_res[i] * t;
    }
}

void ResultReplySub(double* proc_x, double* res, int size, int numprocs, int rank){
    int *RecvNum;
    int *RecvInd;
    int RestLines = size;

    RecvNum = (int*) malloc(numprocs * sizeof(int));
    RecvInd = (int*) malloc(numprocs * sizeof(int));

    RecvInd[0] = 0;
    RecvNum[0] = size / numprocs;
    for (int i = 1; i < numprocs; ++i) {
        RestLines -= RecvNum[i-1];
        RecvNum[i] = RestLines / (numprocs - i);
        RecvInd[i] = RecvInd[i-1] + RecvNum[i-1];
    }
    MPI_Allgatherv(proc_x, RecvNum[rank], MPI_DOUBLE, res, RecvNum, RecvInd,
                   MPI_DOUBLE, MPI_COMM_WORLD);
    free(RecvNum);
    free(RecvInd);
}

int main(int argc, char *argv[]) {
    double *Matrix = NULL;
    double* x0;
    double* vectorB;
    double* result;
    double start_time = 0.0, end_time;
    int numprocs, rank;
    int RestLines;
    int size = DEFAULT;
    if (argc >= 2){
        size = atoi(argv[1]);
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;
    //for each process
    RestLines = size;
    for (int i = 0; i < rank; ++i) {
        RestLines = RestLines - RestLines / (numprocs - i);
    }
    int rowNum = RestLines / (numprocs - rank);
    double* proc_line;
    double* proc_res;
    double* proc_x;
    //memory allocation
    x0 = create_x0(size);
    result = create_x0(size);
    proc_res = (double*) malloc(rowNum * sizeof(double));
    proc_x = (double*) malloc(rowNum * sizeof(double));
    proc_line = (double*) malloc(rowNum * size * sizeof(double));
    vectorB = create_vectorB(size);
    double normB;

    if (rank == ROOT){
        if (size < numprocs) {
            printf("Size of the objects must be greater than number of processes!\n");
            exit(-1);
        }
        Matrix = create_matrix(size);
        normB = norm(vectorB, size);
        for (int i = 1; i < numprocs; ++i) {
            MPI_Send(&normB, 1, MPI_DOUBLE, i, 123, MPI_COMM_WORLD);
        }
        start_time = MPI_Wtime();
    } else {
        MPI_Recv(&normB, 1, MPI_DOUBLE, ROOT, 123, MPI_COMM_WORLD, &status);
    }

    int crit = 0;
    do {
        DataDistribution(Matrix, x0, proc_line, size, rowNum, numprocs, rank);
        ParallelCalculate(proc_line, x0, proc_res, vectorB, size, rowNum, numprocs, rank);
        ResultReply(proc_res, result, size, numprocs, rank);
        if (rank == ROOT){
            crit = criterion(result, normB, size);
            for (int i = 1; i < numprocs; ++i) {
                MPI_Send(&crit, 1, MPI_DOUBLE, i, 123, MPI_COMM_WORLD);
            }
        } else {
            MPI_Recv(&crit, 1, MPI_DOUBLE, ROOT, 123, MPI_COMM_WORLD, &status);
        }
        DataDistributionSub(result, x0, proc_x, proc_res, size, rowNum, numprocs, rank);
        ParallelSubtract(proc_x, result, size, rowNum, numprocs, rank);
        ResultReplySub(proc_x, x0, size, numprocs, rank);
    } while(crit == 0);

    if (rank == ROOT) {
        print_vector(x0, size);
        end_time = MPI_Wtime();
        printf("TIME: %lf\n", end_time - start_time);
    }
    MPI_Finalize();
    FreeProcess(Matrix, vectorB, x0, result, proc_res, proc_line, proc_x, rank);
}