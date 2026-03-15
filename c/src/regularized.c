#include "regularized.h"

#include <complex.h>
#include <math.h>

#include "constants.h"
#include "digamma.h"
#include "fermi_bose.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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

    for (int i = 0; i < n_eps1; ++i) {
        for (int j = 0; j < n_eps2; ++j) {
            double eps1 = epsilon1[i];
            double eps2 = epsilon2[j];

            double _Complex a1 = 0.5 + I * beta * (E2 - eps1) / (2.0 * M_PI);
            double _Complex a2 = 0.5 - I * beta * (E2 - eps2) / (2.0 * M_PI);
            double _Complex a3 = 0.5 + I * beta * (E1 - eps1) / (2.0 * M_PI);
            double _Complex a4 = 0.5 - I * beta * (E1 - eps2) / (2.0 * M_PI);

            double val =
                bose_val / (eps1 - eps2) *
                creal(digamma_c(a1) - digamma_c(a2) - digamma_c(a3) + digamma_c(a4));

            if (!isfinite(val)) {
                val = 0.0;
            }
            out[i * n_eps2 + j] = val;
        }
    }
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

    for (int i = 0; i < n_eps; ++i) {
        double eps = epsilon[i];
        for (int j = 0; j < n_E2; ++j) {
            double e2 = E2[j];
            double _Complex a1 = 0.5 + I * beta * (e2 - eps) / (2.0 * M_PI);
            double _Complex a2 = 0.5 + I * beta * (E1 - eps) / (2.0 * M_PI);

            double val =
                beta / (2.0 * M_PI) * bose_fcn(e2 - E1, T) * cimag(trigamma_c(a1) - trigamma_c(a2));

            if (!isfinite(val)) {
                val = 0.0;
            }
            out[i * n_E2 + j] = val;
        }
    }
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

    for (int i = 0; i < eps_cols; ++i) {
        for (int j = 0; j < n_E2; ++j) {
            if (j >= eps_rows) {
                out[i * n_E2 + j] = 0.0;
                continue;
            }

            double eps = epsilon_matrix[j * eps_cols + i];
            double e2 = E2[j];
            double _Complex a1 = 0.5 + I * beta * (e2 - eps) / (2.0 * M_PI);
            double _Complex a2 = 0.5 + I * beta * (E1 - eps) / (2.0 * M_PI);

            double val =
                beta / (2.0 * M_PI) * bose_fcn(e2 - E1, T) * cimag(trigamma_c(a1) - trigamma_c(a2));

            if (!isfinite(val)) {
                val = 0.0;
            }
            out[i * n_E2 + j] = val;
        }
    }
}
