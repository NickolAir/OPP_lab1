void subtract (double *x, const double *y, int N){
    for (int i = 0; i < N; ++i) {
        x[i] = x[i] - y[i];
    }
}