#ifndef FC_MATRIX_H
#define FC_MATRIX_H

#include "rate.h"

double peq_(int q, double vmode, double T);
double sigma_W(const RateStore *store, int N, int n1, int n2, int q1);
void generate_matrix_W(const RateStore *store, int N,
                       double vmode, double T, double tau,
                       double *W);

#endif
