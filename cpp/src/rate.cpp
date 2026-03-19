#include "rate.hpp"

#include <cmath>
#include <vector>

#include "constants.hpp"
#include "cotunneling.hpp"
#include "fermi_bose.hpp"

namespace fc {

static constexpr double PI = 3.14159265358979323846;

static inline int lead_to_idx(int lead) {
    return (lead == 1) ? 0 : 1;
}

static inline int rate_idx(int N, int n1, int n2, int q1, int lead_idx, int q2) {
    return ((((n1 * 2 + n2) * N + q1) * 2 + lead_idx) * N + q2);
}

static inline double sanitize(double x) {
    return std::isfinite(x) ? x : 0.0;
}

RateStore::RateStore(int N) : N(N), data(4 * N * 2 * N, 0.0) {}

double spin_degeneracy(int n1, int n2) {
    if (n1 == 0) {
        return 2.0;
    }
    return (n2 == 0) ? 1.0 : 2.0;
}

void m_rateW(int n1, int n2, int q1, const int* q2_vec, int nq2,
             double vmode, double alphaL, double alphaR, double lambda,
             double Vsd, double T, double eta, int lead, double Vg,
             const FCCache& fc, double* out) {
    (void)lambda;
    if (nq2 <= 0) {
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
            const double f = fermi(epsilond - static_cast<double>(q2 - q1) * vmode, mu, T);
            const double val = s * gamma / HBAR_EV * fc_val * fc_val * (1.0 - f);
            out[i] = sanitize(val);
        }
    } else if (n1 == 0 && n2 == 1) {
        for (int i = 0; i < nq2; ++i) {
            const int q2 = q2_vec[i];
            const double fc_val = fc_cache_get(fc, q1, q2);
            const double f = fermi(epsilond + static_cast<double>(q2 - q1) * vmode, mu, T);
            const double val = s * gamma / HBAR_EV * fc_val * fc_val * f;
            out[i] = sanitize(val);
        }
    } else if (n1 == 0 && n2 == 0) {
        std::vector<double> single_sum(nq2);
        std::vector<double> double_sum(nq2);

        if (lead == 1) {
            sumMMr(q1, q2_vec, nq2, lambda, muL, muR, vmode, epsilond, T, fc, single_sum.data());
            sumMMMMrs(q1, q2_vec, nq2, lambda, muL, muR, vmode, epsilond, T, fc, double_sum.data());
        } else {
            sumMMr(q1, q2_vec, nq2, lambda, muR, muL, vmode, epsilond, T, fc, single_sum.data());
            sumMMMMrs(q1, q2_vec, nq2, lambda, muR, muL, vmode, epsilond, T, fc, double_sum.data());
        }

        const double prefac = s / (2.0 * PI * HBAR_EV) * gammaL * gammaR;
        for (int i = 0; i < nq2; ++i) {
            out[i] = sanitize(prefac * (single_sum[i] + double_sum[i]));
        }
    } else {
        std::vector<double> single_sum(nq2);
        std::vector<double> double_sum(nq2);

        if (lead == 1) {
            sumMMr11(q1, q2_vec, nq2, lambda, muL, muR, vmode, epsilond, T, fc, single_sum.data());
            sumMMMMrs11(q1, q2_vec, nq2, lambda, muL, muR, vmode, epsilond, T, fc, double_sum.data());
        } else {
            sumMMr11(q1, q2_vec, nq2, lambda, muR, muL, vmode, epsilond, T, fc, single_sum.data());
            sumMMMMrs11(q1, q2_vec, nq2, lambda, muR, muL, vmode, epsilond, T, fc, double_sum.data());
        }

        const double prefac = s / (2.0 * PI * HBAR_EV) * gammaL * gammaR;
        for (int i = 0; i < nq2; ++i) {
            out[i] = sanitize(prefac * (single_sum[i] + double_sum[i]));
        }
    }
}

void calculate_all_rateW(RateStore& store, int N,
                         double vmode, double alphaL, double alphaR,
                         double lambda, double Vsd, double T,
                         double eta, int lead, double Vg, const FCCache& fc) {
    if (N <= 0) {
        return;
    }

    std::vector<int> q2_vec(N);
    std::vector<double> out(N);

    for (int q2 = 0; q2 < N; ++q2) {
        q2_vec[q2] = q2;
    }

    const int lead_idx_val = lead_to_idx(lead);
    for (int n1 = 0; n1 <= 1; ++n1) {
        for (int n2 = 0; n2 <= 1; ++n2) {
            for (int q1 = 0; q1 < N; ++q1) {
                m_rateW(n1, n2, q1, q2_vec.data(), N, vmode, alphaL, alphaR, lambda,
                        Vsd, T, eta, lead, Vg, fc, out.data());
                for (int q2 = 0; q2 < N; ++q2) {
                    const int idx = rate_idx(N, n1, n2, q1, lead_idx_val, q2);
                    store.data[idx] = out[q2];
                }
            }
        }
    }
}

double rateW_from_store(const RateStore& store, int n1, int n2, int q1, int q2, int lead) {
    const int N = store.N;
    if (n1 < 0 || n1 > 1 || n2 < 0 || n2 > 1 || q1 < 0 || q1 >= N || q2 < 0 || q2 >= N) {
        return 0.0;
    }
    const int lid = lead_to_idx(lead);
    const int idx = rate_idx(N, n1, n2, q1, lid, q2);
    return store.data[idx];
}

double rateW_lead_from_store(const RateStore& store, int n1, int n2, int q1, int q2) {
    return rateW_from_store(store, n1, n2, q1, q2, 1) +
           rateW_from_store(store, n1, n2, q1, q2, -1);
}

}
