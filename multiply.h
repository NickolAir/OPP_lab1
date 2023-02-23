#ifndef OPP_LAB1_MULTIPLY_H
#define OPP_LAB1_MULTIPLY_H

double *multiply (double **A, int N, double *x);
void mult_on_vector(double *x, double *y, int N, int rank, int numprocs, double *res);
void mult_on_const(double *x, int N, double t);

#endif //OPP_LAB1_MULTIPLY_H