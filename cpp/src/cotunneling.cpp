#include "cotunneling.hpp"

#include <cmath>
#include <complex>
#include <vector>

#include "regularized.hpp"

namespace fc {

static inline double sanitize(double x) {
    return std::isfinite(x) ? x : 0.0;
}

static double rel_diff_log10(double a, double b) {
    std::complex<double> z(b / a, 0.0);
    if (a == 0.0 && b == 0.0) return 0.0;
    z = std::complex<double>(b, 0.0) / std::complex<double>(a, 0.0);
    return std::abs(std::log(z) / std::log(10.0));
}

static double rel_diff_log(double a, double b) {
    std::complex<double> z = std::complex<double>(b, 0.0) / std::complex<double>(a, 0.0);
    return std::abs(std::log(z));
}

static bool any_greater_than(const double* arr, int n, double threshold) {
    for (int i = 0; i < n; ++i) {
        if (arr[i] > threshold) {
            return true;
        }
    }
    return false;
}

void m_sumMMr(int N, int q1, const int* q2_vec, int nq2,
              double lambda, double muL, double muR,
              double vmode, double epsilond, double T,
              FCCache& fc, double* out) {
    (void)lambda;
    if (N <= 0 || nq2 <= 0) {
        return;
    }

    std::vector<double> MM_sq(N * nq2);
    std::vector<double> Jr(N * nq2);
    std::vector<double> E2_vec(nq2);
    std::vector<double> eps_vec(N);

    for (int j = 0; j < nq2; ++j) {
        E2_vec[j] = muR - (q1 - q2_vec[j]) * vmode;
    }
    for (int r = 0; r < N; ++r) {
        eps_vec[r] = epsilond - (q1 - r) * vmode;
    }

    for (int r = 0; r < N; ++r) {
        double fc_q1r = fc_cache_get(fc, q1, r);
        for (int j = 0; j < nq2; ++j) {
            double fc_q2r = fc_cache_get(fc, q2_vec[j], r);
            double prod = fc_q2r * fc_q1r;
            MM_sq[r * nq2 + j] = prod * prod;
        }
    }

    regularized_J(muL, E2_vec.data(), nq2, eps_vec.data(), N, T, Jr.data());

    for (int j = 0; j < nq2; ++j) {
        double sum = 0.0;
        for (int r = 0; r < N; ++r) {
            double term = MM_sq[r * nq2 + j] * Jr[r * nq2 + j];
            sum += sanitize(term);
        }
        out[j] = sanitize(sum);
    }
}

void m_sumMMr11(int N, int q1, const int* q2_vec, int nq2,
                double lambda, double muL, double muR,
                double vmode, double epsilond, double T,
                FCCache& fc, double* out) {
    (void)lambda;
    if (N <= 0 || nq2 <= 0) {
        return;
    }

    std::vector<double> MM_sq(N * nq2);
    std::vector<double> Jr(N * nq2);
    std::vector<double> E2_vec(nq2);
    std::vector<double> eps_matrix(nq2 * N);

    for (int j = 0; j < nq2; ++j) {
        E2_vec[j] = muR - (q1 - q2_vec[j]) * vmode;
    }

    for (int r = 0; r < N; ++r) {
        double fc_q1r = fc_cache_get(fc, q1, r);
        for (int j = 0; j < nq2; ++j) {
            double fc_q2r = fc_cache_get(fc, q2_vec[j], r);
            double prod = fc_q2r * fc_q1r;
            MM_sq[r * nq2 + j] = prod * prod;
            eps_matrix[j * N + r] = epsilond + (q2_vec[j] - r) * vmode;
        }
    }

    regularized_J_matrix(muL, E2_vec.data(), nq2, eps_matrix.data(), nq2, N, T, Jr.data());

    for (int j = 0; j < nq2; ++j) {
        double sum = 0.0;
        for (int r = 0; r < N; ++r) {
            double term = MM_sq[r * nq2 + j] * Jr[r * nq2 + j];
            sum += sanitize(term);
        }
        out[j] = sanitize(sum);
    }
}

double m_sumMMMMrs(int N, int q1, int q2,
                   double lambda, double muL, double muR,
                   double vmode, double epsilond, double T,
                   FCCache& fc) {
    (void)lambda;
    if (N <= 0) {
        return 0.0;
    }

    std::vector<double> Irs(N * N);
    std::vector<double> eps1(N);
    std::vector<double> eps2(N);

    for (int r = 0; r < N; ++r) {
        eps1[r] = epsilond - (q1 - r) * vmode;
        eps2[r] = epsilond - (q1 - r) * vmode;
    }

    double E2 = muR - (q1 - q2) * vmode;
    regularized_I(muL, E2, eps1.data(), N, eps2.data(), N, T, Irs.data());

    double sum = 0.0;
    for (int r = 0; r < N; ++r) {
        double fc_q2r = fc_cache_get(fc, q2, r);
        double fc_q1r = fc_cache_get(fc, q1, r);
        for (int s = 0; s < N; ++s) {
            if (r == s) {
                continue;
            }
            double fc_q2s = fc_cache_get(fc, q2, s);
            double fc_q1s = fc_cache_get(fc, q1, s);
            double mmmm = fc_q2r * fc_q1r * fc_q2s * fc_q1s;
            double term = mmmm * Irs[r * N + s];
            sum += sanitize(term);
        }
    }

    return sanitize(sum);
}

double m_sumMMMMrs11(int N, int q1, int q2,
                     double lambda, double muL, double muR,
                     double vmode, double epsilond, double T,
                     FCCache& fc) {
    (void)lambda;
    if (N <= 0) {
        return 0.0;
    }

    std::vector<double> Irs(N * N);
    std::vector<double> eps1(N);
    std::vector<double> eps2(N);

    for (int r = 0; r < N; ++r) {
        eps1[r] = epsilond + (q2 - r) * vmode;
        eps2[r] = epsilond + (q2 - r) * vmode;
    }

    double E2 = muR - (q1 - q2) * vmode;
    regularized_I(muL, E2, eps1.data(), N, eps2.data(), N, T, Irs.data());

    double sum = 0.0;
    for (int r = 0; r < N; ++r) {
        double fc_q2r = fc_cache_get(fc, q2, r);
        double fc_q1r = fc_cache_get(fc, q1, r);
        for (int s = 0; s < N; ++s) {
            if (r == s) {
                continue;
            }
            double fc_q2s = fc_cache_get(fc, q2, s);
            double fc_q1s = fc_cache_get(fc, q1, s);
            double mmmm = fc_q2r * fc_q1r * fc_q2s * fc_q1s;
            double term = mmmm * Irs[r * N + s];
            sum += sanitize(term);
        }
    }

    return sanitize(sum);
}

void sumMMr(int q1, const int* q2_vec, int nq2,
            double lambda, double muL, double muR,
            double vmode, double epsilond, double T,
            FCCache& fc, double* out) {
    if (nq2 <= 0) {
        return;
    }

    int tempN = static_cast<int>(std::llround(std::pow(lambda, 2.2) * 3.0));
    constexpr double epsilon_conv = 1e-14;

    std::vector<double> temp_tol(nq2);
    std::vector<double> temp_tol2(nq2);
    std::vector<double> relative_diff(nq2);
    std::vector<int> loc_indices(nq2);

    m_sumMMr(tempN, q1, q2_vec, nq2, lambda, muL, muR, vmode, epsilond, T, fc, temp_tol.data());
    tempN += 5;
    m_sumMMr(tempN, q1, q2_vec, nq2, lambda, muL, muR, vmode, epsilond, T, fc, temp_tol2.data());

    for (int j = 0; j < nq2; ++j) {
        relative_diff[j] = rel_diff_log10(temp_tol[j], temp_tol2[j]);
    }

    while (any_greater_than(relative_diff.data(), nq2, epsilon_conv)) {
        int nloc = 0;
        for (int j = 0; j < nq2; ++j) {
            if (relative_diff[j] > epsilon_conv) {
                loc_indices[nloc++] = j;
            }
        }
        if (nloc == 0) {
            break;
        }

        int step = static_cast<int>(std::llround(tempN * 0.5));
        if (step < 10) step = 10;
        if (step > 20) step = 20;
        tempN += step;

        for (int j = 0; j < nq2; ++j) {
            temp_tol[j] = temp_tol2[j];
        }

        std::vector<int> q2_subset(nloc);
        std::vector<double> partial(nloc);

        for (int k = 0; k < nloc; ++k) {
            q2_subset[k] = q2_vec[loc_indices[k]];
        }

        m_sumMMr(tempN, q1, q2_subset.data(), nloc, lambda, muL, muR, vmode, epsilond, T, fc, partial.data());

        for (int k = 0; k < nloc; ++k) {
            temp_tol2[loc_indices[k]] = partial[k];
        }

        for (int j = 0; j < nq2; ++j) {
            relative_diff[j] = rel_diff_log10(temp_tol[j], temp_tol2[j]);
        }
    }

    for (int j = 0; j < nq2; ++j) {
        out[j] = sanitize(temp_tol2[j]);
    }
}

void sumMMr11(int q1, const int* q2_vec, int nq2,
              double lambda, double muL, double muR,
              double vmode, double epsilond, double T,
              FCCache& fc, double* out) {
    if (nq2 <= 0) {
        return;
    }

    int tempN = static_cast<int>(std::llround(std::pow(lambda, 2.2) * 3.0));
    constexpr double epsilon_conv = 1e-14;

    std::vector<double> temp_tol(nq2);
    std::vector<double> temp_tol2(nq2);
    std::vector<double> relative_diff(nq2);
    std::vector<int> loc_indices(nq2);

    m_sumMMr11(tempN, q1, q2_vec, nq2, lambda, muL, muR, vmode, epsilond, T, fc, temp_tol.data());
    tempN += 5;
    m_sumMMr11(tempN, q1, q2_vec, nq2, lambda, muL, muR, vmode, epsilond, T, fc, temp_tol2.data());

    for (int j = 0; j < nq2; ++j) {
        relative_diff[j] = rel_diff_log10(temp_tol[j], temp_tol2[j]);
    }

    while (any_greater_than(relative_diff.data(), nq2, epsilon_conv)) {
        int nloc = 0;
        for (int j = 0; j < nq2; ++j) {
            if (relative_diff[j] > epsilon_conv) {
                loc_indices[nloc++] = j;
            }
        }
        if (nloc == 0) {
            break;
        }

        int step = static_cast<int>(std::llround(tempN * 0.5));
        if (step < 10) step = 10;
        if (step > 20) step = 20;
        tempN += step;

        for (int j = 0; j < nq2; ++j) {
            temp_tol[j] = temp_tol2[j];
        }

        std::vector<int> q2_subset(nloc);
        std::vector<double> partial(nloc);

        for (int k = 0; k < nloc; ++k) {
            q2_subset[k] = q2_vec[loc_indices[k]];
        }

        m_sumMMr11(tempN, q1, q2_subset.data(), nloc, lambda, muL, muR, vmode, epsilond, T, fc, partial.data());

        for (int k = 0; k < nloc; ++k) {
            temp_tol2[loc_indices[k]] = partial[k];
        }

        for (int j = 0; j < nq2; ++j) {
            relative_diff[j] = rel_diff_log10(temp_tol[j], temp_tol2[j]);
        }
    }

    for (int j = 0; j < nq2; ++j) {
        out[j] = sanitize(temp_tol2[j]);
    }
}

void sumMMMMrs(int q1, const int* q2_vec, int nq2,
               double lambda, double muL, double muR,
               double vmode, double epsilond, double T,
               FCCache& fc, double* out) {
    if (nq2 <= 0) {
        return;
    }

    for (int idx = 0; idx < nq2; ++idx) {
        int q2 = q2_vec[idx];
        int tempN = static_cast<int>(std::llround(lambda * lambda * 4.0));
        constexpr double epsilon_conv = 1e-14;

        double temp_tol = m_sumMMMMrs(tempN, q1, q2, lambda, muL, muR, vmode, epsilond, T, fc);
        tempN += 5;
        double temp_tol2 = m_sumMMMMrs(tempN, q1, q2, lambda, muL, muR, vmode, epsilond, T, fc);

        double relative_diff_val = rel_diff_log(temp_tol, temp_tol2);

        while (std::isfinite(relative_diff_val) && relative_diff_val > epsilon_conv) {
            temp_tol = temp_tol2;
            int step = static_cast<int>(std::llround(tempN * 0.5));
            if (step < 20) step = 20;
            if (step > 40) step = 40;
            tempN += step;
            temp_tol2 = m_sumMMMMrs(tempN, q1, q2, lambda, muL, muR, vmode, epsilond, T, fc);
            relative_diff_val = rel_diff_log(temp_tol, temp_tol2);
        }

        out[idx] = sanitize(temp_tol2);
    }
}

void sumMMMMrs11(int q1, const int* q2_vec, int nq2,
                 double lambda, double muL, double muR,
                 double vmode, double epsilond, double T,
                 FCCache& fc, double* out) {
    if (nq2 <= 0) {
        return;
    }

    for (int idx = 0; idx < nq2; ++idx) {
        int q2 = q2_vec[idx];
        int tempN = static_cast<int>(std::llround(lambda * lambda * 4.0));
        constexpr double epsilon_conv = 1e-14;

        double temp_tol = m_sumMMMMrs11(tempN, q1, q2, lambda, muL, muR, vmode, epsilond, T, fc);
        tempN += 5;
        int secondN = tempN - 1;
        if (secondN < 2) secondN = 2;
        double temp_tol2 = m_sumMMMMrs11(secondN, q1, q2, lambda, muL, muR, vmode, epsilond, T, fc);

        double relative_diff_val = rel_diff_log(temp_tol, temp_tol2);

        while (std::isfinite(relative_diff_val) && relative_diff_val > epsilon_conv) {
            int step = static_cast<int>(std::llround(tempN * 0.5));
            if (step < 20) step = 20;
            if (step > 40) step = 40;
            tempN += step;
            temp_tol = temp_tol2;
            temp_tol2 = m_sumMMMMrs11(tempN, q1, q2, lambda, muL, muR, vmode, epsilond, T, fc);
            relative_diff_val = rel_diff_log(temp_tol, temp_tol2);
        }

        out[idx] = sanitize(temp_tol2);
    }
}

}
