#pragma once

#include <vector>
#include "fc_matrix.hpp"

namespace fc {

struct RateStore {
    int N;
    std::vector<double> data;

    explicit RateStore(int N);
};

double spin_degeneracy(int n1, int n2);

void m_rateW(int n1, int n2, int q1, const int* q2_vec, int nq2,
             double vmode, double alphaL, double alphaR, double lambda,
             double Vsd, double T, double eta, int lead, double Vg,
             const FCCache& fc, double* out);

void calculate_all_rateW(RateStore& store, int N,
                         double vmode, double alphaL, double alphaR,
                         double lambda, double Vsd, double T,
                         double eta, int lead, double Vg, const FCCache& fc);

double rateW_from_store(const RateStore& store, int n1, int n2, int q1, int q2, int lead);
double rateW_lead_from_store(const RateStore& store, int n1, int n2, int q1, int q2);

}
