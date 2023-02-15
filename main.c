#include <mpi.h>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>

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

double *simple_iteration (double **A, double *b, int N, double *x){
    double *new_x = (double*) malloc(sizeof(double) * N);
    for (int i = 0; i < N; ++i) {
        double tmp_sum = x[i];
        for (int j = 0; j < N; ++j) {
                tmp_sum -= A[i][j] * x[j] * t;
        }
        tmp_sum += b[i] * t;
        new_x[i] = tmp_sum;
    }
    return new_x;
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
int criterion (double **A, double *b, int N, double *x_n){
    double *Ax = multiplication(A, N, x_n);
    for (int i = 0; i < N; ++i) {
        Ax[i] -= b[i];
    }
    double result = norm(Ax, N) / norm(b, N);
    if (result < E){
        free(Ax);
        return 1;
    } else {
        free(Ax);
        return 0;
    }
}

void solution (double **Matrix, double *b, double *x, int N){

}

int main(int argc,char *argv[]){
    int erc, numprocs, rank;
    double start_time = 0.0;
    double end_time;

    int N = 3;
    double **Matrix = create_matrix(N);
    double *vector = create_vector(N);
    double *x0 = create_x0(N);

    while (criterion(Matrix, vector, N, x0) == 0){
        x0 = simple_iteration(Matrix, vector, N, x0);
        print_vector(x0, N);
    }

    free(vector);
    free(x0);
    delete_matrix(Matrix, N);

    if (erc = MPI_Init(&argc, &argv)){
        printf("Ошибка запуска, выполнение остановлено\n");
        MPI_Abort(MPI_COMM_WORLD, erc);
    }
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs); // Получение числа процессов
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Получение номера процесса
    int range = N / numprocs;
    if (rank == 0){
        start_time = MPI_Wtime();
        if (N % numprocs == 0){
            for (int to_thread = 1; to_thread < numprocs; to_thread++){
                MPI_Send(&range, 1, MPI_INT, to_thread, 0, MPI_COMM_WORLD);
            }
        } else {
            for (int to_thread = 1; to_thread < numprocs - 1; to_thread++){
                MPI_Send(&range, 1, MPI_INT, to_thread, 0, MPI_COMM_WORLD);
            }
            range += N - range * numprocs;
            MPI_Send(&range, 1, MPI_INT, numprocs - 1, 0, MPI_COMM_WORLD);
        }

        start_time = MPI_Wtime();// Начинаем считать время выполнения
    }
    //status - структура, определенная в MPI которая хранит информацию о пересылке и статус ее завершения.
    MPI_Status status;

    //code
    if (rank == 0){
        end_time = MPI_Wtime();
        printf("%lf\n", (end_time-start_time)*1000);
    }
    MPI_Finalize(); // Завершение работы MPI
    return 0;
}