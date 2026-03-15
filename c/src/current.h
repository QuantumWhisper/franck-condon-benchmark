#ifndef FC_CURRENT_H
#define FC_CURRENT_H

#include "rate.h"

typedef struct {
    double I_tol;
    double I_seq;
    double I_cot;
} CurrentResult;

CurrentResult current_from_rate_equations(
    double Vsd, int N, double vmode, double alphaL, double alphaR,
    double lambda, double T, double eta, double Vg, double tau,
    const RateStore *store);

#endif
