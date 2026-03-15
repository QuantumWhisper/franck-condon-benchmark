#include "matrix.h"

#include <math.h>
#include <stddef.h>

#include "constants.h"

double peq_(int q, double vmode, double T) {
    const double beta = 1.0 / (KB_EV * T);
    return exp(-(double)q * vmode * beta) * (1.0 - exp(-vmode * beta));
}

double sigma_W(const RateStore *store, int N, int n1, int n2, int q1) {
    double s = 0.0;
    for (int q2 = 0; q2 < N; ++q2) {
        s += rateW_lead_from_store(store, n1, n2, q1, q2);
    }
    return s;
}

void generate_matrix_W(const RateStore *store, int N,
                       double vmode, double T, double tau,
                       double *W) {
    if (store == NULL || W == NULL || N <= 0) {
        return;
    }

    const int dim = 2 * N;
    const double inv_tau = 1.0 / tau;

    for (int ii = 0; ii < N; ++ii) {
        for (int jj = 0; jj < dim; ++jj) {
            const int modjj = jj % N;
            const int p_q_index = modjj;
            double val = 0.0;

            if (jj < N) {
                if (ii == modjj) {
                    val = rateW_lead_from_store(store, 0, 0, p_q_index, p_q_index) -
                          sigma_W(store, N, 0, 0, p_q_index) -
                          sigma_W(store, N, 0, 1, p_q_index) -
                          inv_tau + peq_(p_q_index, vmode, T) * inv_tau;
                } else {
                    val = rateW_lead_from_store(store, 0, 0, p_q_index, ii) +
                          peq_(ii, vmode, T) * inv_tau;
                }
            } else {
                val = rateW_lead_from_store(store, 1, 0, p_q_index, ii);
            }

            W[ii * dim + jj] = val;
        }
    }

    for (int ii = 0; ii < N; ++ii) {
        const int row = N + ii;
        for (int jj = 0; jj < dim; ++jj) {
            const int modjj = jj % N;
            double val = 0.0;

            if (jj < N) {
                val = rateW_lead_from_store(store, 0, 1, modjj, ii);
            } else {
                if (modjj == ii) {
                    val = rateW_lead_from_store(store, 1, 1, modjj, ii) -
                          sigma_W(store, N, 1, 0, ii) -
                          sigma_W(store, N, 1, 1, ii) -
                          inv_tau + peq_(ii, vmode, T) * inv_tau;
                } else {
                    val = rateW_lead_from_store(store, 1, 1, modjj, ii) +
                          peq_(ii, vmode, T) * inv_tau;
                }
            }

            W[row * dim + jj] = val;
        }
    }
}
