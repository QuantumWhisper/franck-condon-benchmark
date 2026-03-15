#include "regularized.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>

#include "constants.h"
#include "digamma.h"
#include "fermi_bose.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
 * Key optimization: digamma arguments in regularized_I factor into
 * row-only and column-only terms:
 *   a1[i] = 0.5 + I*beta*(E2-eps1[i])/(2pi)  -- row only
 *   a3[i] = 0.5 + I*beta*(E1-eps1[i])/(2pi)  -- row only
 *   a2[j] = 0.5 - I*beta*(E2-eps2[j])/(2pi)  -- col only
 *   a4[j] = 0.5 - I*beta*(E1-eps2[j])/(2pi)  -- col only
 *
 * Precompute digamma per row/col: 4*N calls instead of 4*N*N.
 */
void regularized_I(double E1,
                   double E2,
                   const double *epsilon1,
                   int n_eps1,
                   const double *epsilon2,
                   int n_eps2,
                   double T,
                   double *out) {
    if (epsilon1 == 0 || epsilon2 == 0 || out == 0 || n_eps1 <= 0 || n_eps2 <= 0) {
        return;
    }

    double beta = 1.0 / (KB_EV * T);
    double bose_val = bose_fcn(E2 - E1, T);
    double coeff = beta / (2.0 * M_PI);

    double _Complex *row_diff = (double _Complex *)malloc((size_t)n_eps1 * sizeof(double _Complex));
    double _Complex *col_diff = (double _Complex *)malloc((size_t)n_eps2 * sizeof(double _Complex));

    for (int i = 0; i < n_eps1; ++i) {
        double _Complex a1 = 0.5 + I * coeff * (E2 - epsilon1[i]);
        double _Complex a3 = 0.5 + I * coeff * (E1 - epsilon1[i]);
        row_diff[i] = digamma_c(a1) - digamma_c(a3);
    }

    for (int j = 0; j < n_eps2; ++j) {
        double _Complex a2 = 0.5 - I * coeff * (E2 - epsilon2[j]);
        double _Complex a4 = 0.5 - I * coeff * (E1 - epsilon2[j]);
        col_diff[j] = digamma_c(a2) - digamma_c(a4);
    }

    for (int i = 0; i < n_eps1; ++i) {
        for (int j = 0; j < n_eps2; ++j) {
            double denom = epsilon1[i] - epsilon2[j];
            double val = bose_val / denom * creal(row_diff[i] - col_diff[j]);
            out[i * n_eps2 + j] = isfinite(val) ? val : 0.0;
        }
    }

    free(row_diff);
    free(col_diff);
}

void regularized_J(double E1,
                   const double *E2,
                   int n_E2,
                   const double *epsilon,
                   int n_eps,
                   double T,
                   double *out) {
    if (E2 == 0 || epsilon == 0 || out == 0 || n_E2 <= 0 || n_eps <= 0) {
        return;
    }

    double beta = 1.0 / (KB_EV * T);
    double coeff = beta / (2.0 * M_PI);

    double *bose_vals = (double *)malloc((size_t)n_E2 * sizeof(double));
    double _Complex *trig_a2 = (double _Complex *)malloc((size_t)n_eps * sizeof(double _Complex));

    for (int j = 0; j < n_E2; ++j) {
        bose_vals[j] = bose_fcn(E2[j] - E1, T);
    }

    for (int i = 0; i < n_eps; ++i) {
        double _Complex a2 = 0.5 + I * coeff * (E1 - epsilon[i]);
        trig_a2[i] = trigamma_c(a2);
    }

    for (int i = 0; i < n_eps; ++i) {
        for (int j = 0; j < n_E2; ++j) {
            double _Complex a1 = 0.5 + I * coeff * (E2[j] - epsilon[i]);
            double val = coeff * bose_vals[j] * cimag(trigamma_c(a1) - trig_a2[i]);
            out[i * n_E2 + j] = isfinite(val) ? val : 0.0;
        }
    }

    free(bose_vals);
    free(trig_a2);
}

void regularized_J_matrix(double E1,
                          const double *E2,
                          int n_E2,
                          const double *epsilon_matrix,
                          int eps_rows,
                          int eps_cols,
                          double T,
                          double *out) {
    if (E2 == 0 || epsilon_matrix == 0 || out == 0 || n_E2 <= 0 || eps_rows <= 0 || eps_cols <= 0) {
        return;
    }

    double beta = 1.0 / (KB_EV * T);
    double coeff = beta / (2.0 * M_PI);

    double *bose_vals = (double *)malloc((size_t)n_E2 * sizeof(double));
    for (int j = 0; j < n_E2; ++j) {
        bose_vals[j] = bose_fcn(E2[j] - E1, T);
    }

    for (int i = 0; i < eps_cols; ++i) {
        for (int j = 0; j < n_E2; ++j) {
            if (j >= eps_rows) {
                out[i * n_E2 + j] = 0.0;
                continue;
            }

            double eps = epsilon_matrix[j * eps_cols + i];
            double _Complex a1 = 0.5 + I * coeff * (E2[j] - eps);
            double _Complex a2 = 0.5 + I * coeff * (E1 - eps);

            double val = coeff * bose_vals[j] * cimag(trigamma_c(a1) - trigamma_c(a2));
            out[i * n_E2 + j] = isfinite(val) ? val : 0.0;
        }
    }

    free(bose_vals);
}
