#include "rate.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#include "constants.h"
#include "cotunneling.h"
#include "fermi_bose.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static inline int lead_to_idx(int lead) {
    return (lead == 1) ? 0 : 1;
}

static inline int rate_idx(int N, int n1, int n2, int q1, int lead_idx, int q2) {
    return ((((n1 * 2 + n2) * N + q1) * 2 + lead_idx) * N + q2);
}

static double sanitize(double x) {
    return isfinite(x) ? x : 0.0;
}

void rate_store_init(RateStore *store, int N) {
    if (store == NULL || N <= 0) {
        return;
    }
    store->N = N;
    store->data = (double *)calloc((size_t)(4 * N * 2 * N), sizeof(double));
}

void rate_store_free(RateStore *store) {
    if (store == NULL) {
        return;
    }
    free(store->data);
    store->data = NULL;
    store->N = 0;
}

double spin_degeneracy(int n1, int n2) {
    if (n1 == 0) {
        return 2.0;
    }
    return (n2 == 0) ? 1.0 : 2.0;
}

void m_rateW(int n1, int n2, int q1, const int *q2_vec, int nq2,
             double vmode, double alphaL, double alphaR, double lambda,
             double Vsd, double T, double eta, int lead, double Vg,
             FCCache *fc, double *out) {
    (void)lambda;
    if (q2_vec == NULL || out == NULL || fc == NULL || nq2 <= 0) {
        return;
    }

    const double gammaL = alphaL * vmode;
    const double gammaR = alphaR * vmode;
    const double epsilond = 0.0 + Vg;
    const double muL = eta * Vsd;
    const double muR = -(1.0 - eta) * Vsd;

    const double s = spin_degeneracy(n1, n2);
    const double gamma = (lead == 1) ? gammaL : gammaR;
    const double mu = (lead == 1) ? muL : muR;

    if (lead != 1 && lead != -1) {
        for (int i = 0; i < nq2; ++i) {
            out[i] = 0.0;
        }
        return;
    }

    if (n1 == 1 && n2 == 0) {
        for (int i = 0; i < nq2; ++i) {
            const int q2 = q2_vec[i];
            const double fc_val = fc_cache_get(fc, q1, q2);
            const double f = fermi(epsilond - (double)(q2 - q1) * vmode, mu, T);
            const double val = s * gamma / HBAR_EV * fc_val * fc_val * (1.0 - f);
            out[i] = sanitize(val);
        }
    } else if (n1 == 0 && n2 == 1) {
        for (int i = 0; i < nq2; ++i) {
            const int q2 = q2_vec[i];
            const double fc_val = fc_cache_get(fc, q1, q2);
            const double f = fermi(epsilond + (double)(q2 - q1) * vmode, mu, T);
            const double val = s * gamma / HBAR_EV * fc_val * fc_val * f;
            out[i] = sanitize(val);
        }
    } else if (n1 == 0 && n2 == 0) {
        double *single_sum = (double *)malloc((size_t)nq2 * sizeof(double));
        double *double_sum = (double *)malloc((size_t)nq2 * sizeof(double));
        if (single_sum == NULL || double_sum == NULL) {
            free(single_sum);
            free(double_sum);
            for (int i = 0; i < nq2; ++i) {
                out[i] = 0.0;
            }
            return;
        }

        if (lead == 1) {
            sumMMr(q1, q2_vec, nq2, lambda, muL, muR, vmode, epsilond, T, fc, single_sum);
            sumMMMMrs(q1, q2_vec, nq2, lambda, muL, muR, vmode, epsilond, T, fc, double_sum);
        } else {
            sumMMr(q1, q2_vec, nq2, lambda, muR, muL, vmode, epsilond, T, fc, single_sum);
            sumMMMMrs(q1, q2_vec, nq2, lambda, muR, muL, vmode, epsilond, T, fc, double_sum);
        }

        const double prefac = s / (2.0 * M_PI * HBAR_EV) * gammaL * gammaR;
        for (int i = 0; i < nq2; ++i) {
            out[i] = sanitize(prefac * (single_sum[i] + double_sum[i]));
        }

        free(single_sum);
        free(double_sum);
    } else {
        double *single_sum = (double *)malloc((size_t)nq2 * sizeof(double));
        double *double_sum = (double *)malloc((size_t)nq2 * sizeof(double));
        if (single_sum == NULL || double_sum == NULL) {
            free(single_sum);
            free(double_sum);
            for (int i = 0; i < nq2; ++i) {
                out[i] = 0.0;
            }
            return;
        }

        if (lead == 1) {
            sumMMr11(q1, q2_vec, nq2, lambda, muL, muR, vmode, epsilond, T, fc, single_sum);
            sumMMMMrs11(q1, q2_vec, nq2, lambda, muL, muR, vmode, epsilond, T, fc, double_sum);
        } else {
            sumMMr11(q1, q2_vec, nq2, lambda, muR, muL, vmode, epsilond, T, fc, single_sum);
            sumMMMMrs11(q1, q2_vec, nq2, lambda, muR, muL, vmode, epsilond, T, fc, double_sum);
        }

        const double prefac = s / (2.0 * M_PI * HBAR_EV) * gammaL * gammaR;
        for (int i = 0; i < nq2; ++i) {
            out[i] = sanitize(prefac * (single_sum[i] + double_sum[i]));
        }

        free(single_sum);
        free(double_sum);
    }
}

void calculate_all_rateW(RateStore *store, int N,
                         double vmode, double alphaL, double alphaR,
                         double lambda, double Vsd, double T,
                         double eta, int lead, double Vg, FCCache *fc) {
    if (store == NULL || store->data == NULL || fc == NULL || N <= 0) {
        return;
    }

    int *q2_vec = (int *)malloc((size_t)N * sizeof(int));
    double *out = (double *)malloc((size_t)N * sizeof(double));
    if (q2_vec == NULL || out == NULL) {
        free(q2_vec);
        free(out);
        return;
    }

    for (int q2 = 0; q2 < N; ++q2) {
        q2_vec[q2] = q2;
    }

    const int lead_idx = lead_to_idx(lead);
    for (int n1 = 0; n1 <= 1; ++n1) {
        for (int n2 = 0; n2 <= 1; ++n2) {
            for (int q1 = 0; q1 < N; ++q1) {
                m_rateW(n1, n2, q1, q2_vec, N, vmode, alphaL, alphaR, lambda,
                        Vsd, T, eta, lead, Vg, fc, out);
                for (int q2 = 0; q2 < N; ++q2) {
                    const int idx = rate_idx(N, n1, n2, q1, lead_idx, q2);
                    store->data[idx] = out[q2];
                }
            }
        }
    }

    free(q2_vec);
    free(out);
}

double rateW_from_store(const RateStore *store, int n1, int n2, int q1, int q2, int lead) {
    if (store == NULL || store->data == NULL || store->N <= 0) {
        return 0.0;
    }
    const int N = store->N;
    if (n1 < 0 || n1 > 1 || n2 < 0 || n2 > 1 || q1 < 0 || q1 >= N || q2 < 0 || q2 >= N) {
        return 0.0;
    }
    const int lead_idx = lead_to_idx(lead);
    const int idx = rate_idx(N, n1, n2, q1, lead_idx, q2);
    return store->data[idx];
}

double rateW_lead_from_store(const RateStore *store, int n1, int n2, int q1, int q2) {
    const double wl = rateW_from_store(store, n1, n2, q1, q2, 1);
    const double wr = rateW_from_store(store, n1, n2, q1, q2, -1);
    return wl + wr;
}
