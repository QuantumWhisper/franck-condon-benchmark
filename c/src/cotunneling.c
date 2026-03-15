#include "cotunneling.h"

#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#include "regularized.h"

static double sanitize(double x) {
    return isfinite(x) ? x : 0.0;
}

static double rel_diff_log10(double a, double b) {
    double _Complex z = b / a;
    return cabs(clog(z) / log(10.0));
}

static double rel_diff_log(double a, double b) {
    double _Complex z = b / a;
    return cabs(clog(z));
}

static int any_greater_than(const double *arr, int n, double threshold) {
    for (int i = 0; i < n; ++i) {
        if (arr[i] > threshold) {
            return 1;
        }
    }
    return 0;
}

void m_sumMMr(int N,
              int q1,
              const int *q2_vec,
              int nq2,
              double lambda,
              double muL,
              double muR,
              double vmode,
              double epsilond,
              double T,
              FCCache *fc,
              double *out) {
    (void)lambda;
    if (q2_vec == 0 || out == 0 || fc == 0 || N <= 0 || nq2 <= 0) {
        return;
    }

    double *MM_sq = (double *)malloc((size_t)N * (size_t)nq2 * sizeof(double));
    double *Jr = (double *)malloc((size_t)N * (size_t)nq2 * sizeof(double));
    double *E2_vec = (double *)malloc((size_t)nq2 * sizeof(double));
    double *eps_vec = (double *)malloc((size_t)N * sizeof(double));

    if (MM_sq == 0 || Jr == 0 || E2_vec == 0 || eps_vec == 0) {
        free(MM_sq);
        free(Jr);
        free(E2_vec);
        free(eps_vec);
        for (int j = 0; j < nq2; ++j) {
            out[j] = 0.0;
        }
        return;
    }

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

    regularized_J(muL, E2_vec, nq2, eps_vec, N, T, Jr);

    for (int j = 0; j < nq2; ++j) {
        double sum = 0.0;
        for (int r = 0; r < N; ++r) {
            double term = MM_sq[r * nq2 + j] * Jr[r * nq2 + j];
            sum += sanitize(term);
        }
        out[j] = sanitize(sum);
    }

    free(MM_sq);
    free(Jr);
    free(E2_vec);
    free(eps_vec);
}

void m_sumMMr11(int N,
                int q1,
                const int *q2_vec,
                int nq2,
                double lambda,
                double muL,
                double muR,
                double vmode,
                double epsilond,
                double T,
                FCCache *fc,
                double *out) {
    (void)lambda;
    if (q2_vec == 0 || out == 0 || fc == 0 || N <= 0 || nq2 <= 0) {
        return;
    }

    double *MM_sq = (double *)malloc((size_t)N * (size_t)nq2 * sizeof(double));
    double *Jr = (double *)malloc((size_t)N * (size_t)nq2 * sizeof(double));
    double *E2_vec = (double *)malloc((size_t)nq2 * sizeof(double));
    double *eps_matrix = (double *)malloc((size_t)nq2 * (size_t)N * sizeof(double));

    if (MM_sq == 0 || Jr == 0 || E2_vec == 0 || eps_matrix == 0) {
        free(MM_sq);
        free(Jr);
        free(E2_vec);
        free(eps_matrix);
        for (int j = 0; j < nq2; ++j) {
            out[j] = 0.0;
        }
        return;
    }

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

    regularized_J_matrix(muL, E2_vec, nq2, eps_matrix, nq2, N, T, Jr);

    for (int j = 0; j < nq2; ++j) {
        double sum = 0.0;
        for (int r = 0; r < N; ++r) {
            double term = MM_sq[r * nq2 + j] * Jr[r * nq2 + j];
            sum += sanitize(term);
        }
        out[j] = sanitize(sum);
    }

    free(MM_sq);
    free(Jr);
    free(E2_vec);
    free(eps_matrix);
}

double m_sumMMMMrs(int N,
                   int q1,
                   int q2,
                   double lambda,
                   double muL,
                   double muR,
                   double vmode,
                   double epsilond,
                   double T,
                   FCCache *fc) {
    (void)lambda;
    if (fc == 0 || N <= 0) {
        return 0.0;
    }

    double *Irs = (double *)malloc((size_t)N * (size_t)N * sizeof(double));
    double *eps1 = (double *)malloc((size_t)N * sizeof(double));
    double *eps2 = (double *)malloc((size_t)N * sizeof(double));

    if (Irs == 0 || eps1 == 0 || eps2 == 0) {
        free(Irs);
        free(eps1);
        free(eps2);
        return 0.0;
    }

    for (int r = 0; r < N; ++r) {
        eps1[r] = epsilond - (q1 - r) * vmode;
        eps2[r] = epsilond - (q1 - r) * vmode;
    }

    double E2 = muR - (q1 - q2) * vmode;
    regularized_I(muL, E2, eps1, N, eps2, N, T, Irs);

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

    free(Irs);
    free(eps1);
    free(eps2);

    return sanitize(sum);
}

double m_sumMMMMrs11(int N,
                     int q1,
                     int q2,
                     double lambda,
                     double muL,
                     double muR,
                     double vmode,
                     double epsilond,
                     double T,
                     FCCache *fc) {
    (void)lambda;
    if (fc == 0 || N <= 0) {
        return 0.0;
    }

    double *Irs = (double *)malloc((size_t)N * (size_t)N * sizeof(double));
    double *eps1 = (double *)malloc((size_t)N * sizeof(double));
    double *eps2 = (double *)malloc((size_t)N * sizeof(double));

    if (Irs == 0 || eps1 == 0 || eps2 == 0) {
        free(Irs);
        free(eps1);
        free(eps2);
        return 0.0;
    }

    for (int r = 0; r < N; ++r) {
        eps1[r] = epsilond + (q2 - r) * vmode;
        eps2[r] = epsilond + (q2 - r) * vmode;
    }

    double E2 = muR - (q1 - q2) * vmode;
    regularized_I(muL, E2, eps1, N, eps2, N, T, Irs);

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

    free(Irs);
    free(eps1);
    free(eps2);

    return sanitize(sum);
}

void sumMMr(int q1,
            const int *q2_vec,
            int nq2,
            double lambda,
            double muL,
            double muR,
            double vmode,
            double epsilond,
            double T,
            FCCache *fc,
            double *out) {
    if (q2_vec == 0 || out == 0 || fc == 0 || nq2 <= 0) {
        return;
    }

    int tempN = (int)llround(pow(lambda, 2.2) * 3.0);
    const double epsilon_conv = 1e-14;

    double *temp_tol = (double *)malloc((size_t)nq2 * sizeof(double));
    double *temp_tol2 = (double *)malloc((size_t)nq2 * sizeof(double));
    double *relative_diff = (double *)malloc((size_t)nq2 * sizeof(double));
    int *loc_indices = (int *)malloc((size_t)nq2 * sizeof(int));

    if (temp_tol == 0 || temp_tol2 == 0 || relative_diff == 0 || loc_indices == 0) {
        free(temp_tol);
        free(temp_tol2);
        free(relative_diff);
        free(loc_indices);
        for (int j = 0; j < nq2; ++j) {
            out[j] = 0.0;
        }
        return;
    }

    m_sumMMr(tempN, q1, q2_vec, nq2, lambda, muL, muR, vmode, epsilond, T, fc, temp_tol);
    tempN += 5;
    m_sumMMr(tempN, q1, q2_vec, nq2, lambda, muL, muR, vmode, epsilond, T, fc, temp_tol2);

    for (int j = 0; j < nq2; ++j) {
        relative_diff[j] = rel_diff_log10(temp_tol[j], temp_tol2[j]);
    }

    while (any_greater_than(relative_diff, nq2, epsilon_conv)) {
        int nloc = 0;
        for (int j = 0; j < nq2; ++j) {
            if (relative_diff[j] > epsilon_conv) {
                loc_indices[nloc++] = j;
            }
        }
        if (nloc == 0) {
            break;
        }

        int step = (int)llround(tempN * 0.5);
        if (step < 10) {
            step = 10;
        }
        if (step > 20) {
            step = 20;
        }
        tempN += step;

        for (int j = 0; j < nq2; ++j) {
            temp_tol[j] = temp_tol2[j];
        }

        int *q2_subset = (int *)malloc((size_t)nloc * sizeof(int));
        double *partial = (double *)malloc((size_t)nloc * sizeof(double));
        if (q2_subset == 0 || partial == 0) {
            free(q2_subset);
            free(partial);
            break;
        }

        for (int k = 0; k < nloc; ++k) {
            q2_subset[k] = q2_vec[loc_indices[k]];
        }

        m_sumMMr(tempN, q1, q2_subset, nloc, lambda, muL, muR, vmode, epsilond, T, fc, partial);

        for (int k = 0; k < nloc; ++k) {
            temp_tol2[loc_indices[k]] = partial[k];
        }

        free(q2_subset);
        free(partial);

        for (int j = 0; j < nq2; ++j) {
            relative_diff[j] = rel_diff_log10(temp_tol[j], temp_tol2[j]);
        }
    }

    for (int j = 0; j < nq2; ++j) {
        out[j] = sanitize(temp_tol2[j]);
    }

    free(temp_tol);
    free(temp_tol2);
    free(relative_diff);
    free(loc_indices);
}

void sumMMr11(int q1,
              const int *q2_vec,
              int nq2,
              double lambda,
              double muL,
              double muR,
              double vmode,
              double epsilond,
              double T,
              FCCache *fc,
              double *out) {
    if (q2_vec == 0 || out == 0 || fc == 0 || nq2 <= 0) {
        return;
    }

    int tempN = (int)llround(pow(lambda, 2.2) * 3.0);
    const double epsilon_conv = 1e-14;

    double *temp_tol = (double *)malloc((size_t)nq2 * sizeof(double));
    double *temp_tol2 = (double *)malloc((size_t)nq2 * sizeof(double));
    double *relative_diff = (double *)malloc((size_t)nq2 * sizeof(double));
    int *loc_indices = (int *)malloc((size_t)nq2 * sizeof(int));

    if (temp_tol == 0 || temp_tol2 == 0 || relative_diff == 0 || loc_indices == 0) {
        free(temp_tol);
        free(temp_tol2);
        free(relative_diff);
        free(loc_indices);
        for (int j = 0; j < nq2; ++j) {
            out[j] = 0.0;
        }
        return;
    }

    m_sumMMr11(tempN, q1, q2_vec, nq2, lambda, muL, muR, vmode, epsilond, T, fc, temp_tol);
    tempN += 5;
    m_sumMMr11(tempN, q1, q2_vec, nq2, lambda, muL, muR, vmode, epsilond, T, fc, temp_tol2);

    for (int j = 0; j < nq2; ++j) {
        relative_diff[j] = rel_diff_log10(temp_tol[j], temp_tol2[j]);
    }

    while (any_greater_than(relative_diff, nq2, epsilon_conv)) {
        int nloc = 0;
        for (int j = 0; j < nq2; ++j) {
            if (relative_diff[j] > epsilon_conv) {
                loc_indices[nloc++] = j;
            }
        }
        if (nloc == 0) {
            break;
        }

        int step = (int)llround(tempN * 0.5);
        if (step < 10) {
            step = 10;
        }
        if (step > 20) {
            step = 20;
        }
        tempN += step;

        for (int j = 0; j < nq2; ++j) {
            temp_tol[j] = temp_tol2[j];
        }

        int *q2_subset = (int *)malloc((size_t)nloc * sizeof(int));
        double *partial = (double *)malloc((size_t)nloc * sizeof(double));
        if (q2_subset == 0 || partial == 0) {
            free(q2_subset);
            free(partial);
            break;
        }

        for (int k = 0; k < nloc; ++k) {
            q2_subset[k] = q2_vec[loc_indices[k]];
        }

        m_sumMMr11(tempN, q1, q2_subset, nloc, lambda, muL, muR, vmode, epsilond, T, fc, partial);

        for (int k = 0; k < nloc; ++k) {
            temp_tol2[loc_indices[k]] = partial[k];
        }

        free(q2_subset);
        free(partial);

        for (int j = 0; j < nq2; ++j) {
            relative_diff[j] = rel_diff_log10(temp_tol[j], temp_tol2[j]);
        }
    }

    for (int j = 0; j < nq2; ++j) {
        out[j] = sanitize(temp_tol2[j]);
    }

    free(temp_tol);
    free(temp_tol2);
    free(relative_diff);
    free(loc_indices);
}

void sumMMMMrs(int q1,
               const int *q2_vec,
               int nq2,
               double lambda,
               double muL,
               double muR,
               double vmode,
               double epsilond,
               double T,
               FCCache *fc,
               double *out) {
    if (q2_vec == 0 || out == 0 || fc == 0 || nq2 <= 0) {
        return;
    }

    for (int idx = 0; idx < nq2; ++idx) {
        int q2 = q2_vec[idx];
        int tempN = (int)llround(lambda * lambda * 4.0);
        const double epsilon_conv = 1e-14;

        double temp_tol = m_sumMMMMrs(tempN, q1, q2, lambda, muL, muR, vmode, epsilond, T, fc);
        tempN += 5;
        double temp_tol2 = m_sumMMMMrs(tempN, q1, q2, lambda, muL, muR, vmode, epsilond, T, fc);

        double relative_diff = rel_diff_log(temp_tol, temp_tol2);

        while (isfinite(relative_diff) && relative_diff > epsilon_conv) {
            temp_tol = temp_tol2;
            int step = (int)llround(tempN * 0.5);
            if (step < 20) {
                step = 20;
            }
            if (step > 40) {
                step = 40;
            }
            tempN += step;
            temp_tol2 = m_sumMMMMrs(tempN, q1, q2, lambda, muL, muR, vmode, epsilond, T, fc);
            relative_diff = rel_diff_log(temp_tol, temp_tol2);
        }

        out[idx] = sanitize(temp_tol2);
    }
}

void sumMMMMrs11(int q1,
                 const int *q2_vec,
                 int nq2,
                 double lambda,
                 double muL,
                 double muR,
                 double vmode,
                 double epsilond,
                 double T,
                 FCCache *fc,
                 double *out) {
    if (q2_vec == 0 || out == 0 || fc == 0 || nq2 <= 0) {
        return;
    }

    for (int idx = 0; idx < nq2; ++idx) {
        int q2 = q2_vec[idx];
        int tempN = (int)llround(lambda * lambda * 4.0);
        const double epsilon_conv = 1e-14;

        double temp_tol = m_sumMMMMrs11(tempN, q1, q2, lambda, muL, muR, vmode, epsilond, T, fc);
        tempN += 5;
        int secondN = tempN - 1;
        if (secondN < 2) {
            secondN = 2;
        }
        double temp_tol2 =
            m_sumMMMMrs11(secondN, q1, q2, lambda, muL, muR, vmode, epsilond, T, fc);

        double relative_diff = rel_diff_log(temp_tol, temp_tol2);

        while (isfinite(relative_diff) && relative_diff > epsilon_conv) {
            int step = (int)llround(tempN * 0.5);
            if (step < 20) {
                step = 20;
            }
            if (step > 40) {
                step = 40;
            }
            tempN += step;
            temp_tol = temp_tol2;
            temp_tol2 = m_sumMMMMrs11(tempN, q1, q2, lambda, muL, muR, vmode, epsilond, T, fc);
            relative_diff = rel_diff_log(temp_tol, temp_tol2);
        }

        out[idx] = sanitize(temp_tol2);
    }
}
