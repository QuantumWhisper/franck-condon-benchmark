#ifndef FC_RATE_H
#define FC_RATE_H

#include "fc_matrix.h"

typedef struct {
    int N;
    double *data;
} RateStore;

void rate_store_init(RateStore *store, int N);
void rate_store_free(RateStore *store);

double spin_degeneracy(int n1, int n2);

void m_rateW(int n1, int n2, int q1, const int *q2_vec, int nq2,
             double vmode, double alphaL, double alphaR, double lambda,
             double Vsd, double T, double eta, int lead, double Vg,
             FCCache *fc, double *out);

void calculate_all_rateW(RateStore *store, int N,
                         double vmode, double alphaL, double alphaR,
                         double lambda, double Vsd, double T,
                         double eta, int lead, double Vg, FCCache *fc);

double rateW_from_store(const RateStore *store, int n1, int n2, int q1, int q2, int lead);
double rateW_lead_from_store(const RateStore *store, int n1, int n2, int q1, int q2);

#endif
