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
    double new_x[N];
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

double *multiplication (double **A, int N, double *x){
    double new_x[N];
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
    double *res = multiplication(Matrix, N, vector);
    print_matrix(Matrix, N);

    print_vector(res, N);

/*    if (rc = MPI_Init(&argc, &argv))
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
    return 0;*/
}