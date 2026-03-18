#include "regularized.hpp"

#include <cmath>
#include <complex>
#include <vector>

#include "constants.hpp"
#include "digamma.hpp"
#include "fermi_bose.hpp"

namespace fc {

static constexpr double PI = 3.14159265358979323846;

void regularized_I(double E1, double E2,
                   const double* epsilon1, int n_eps1,
                   const double* epsilon2, int n_eps2,
                   double T, double* out) {
    if (n_eps1 <= 0 || n_eps2 <= 0) {
        return;
    }

    double beta = 1.0 / (KB_EV * T);
    double bose_val = bose_fcn(E2 - E1, T);
    double coeff = beta / (2.0 * PI);

    std::vector<Complex> row_diff(n_eps1);
    std::vector<Complex> col_diff(n_eps2);

    for (int i = 0; i < n_eps1; ++i) {
        Complex a1 = 0.5 + Complex(0.0, 1.0) * coeff * (E2 - epsilon1[i]);
        Complex a3 = 0.5 + Complex(0.0, 1.0) * coeff * (E1 - epsilon1[i]);
        row_diff[i] = digamma(a1) - digamma(a3);
    }

    for (int j = 0; j < n_eps2; ++j) {
        Complex a2 = 0.5 - Complex(0.0, 1.0) * coeff * (E2 - epsilon2[j]);
        Complex a4 = 0.5 - Complex(0.0, 1.0) * coeff * (E1 - epsilon2[j]);
        col_diff[j] = digamma(a2) - digamma(a4);
    }

    for (int i = 0; i < n_eps1; ++i) {
        for (int j = 0; j < n_eps2; ++j) {
            double denom = epsilon1[i] - epsilon2[j];
            double val = bose_val / denom * (row_diff[i] - col_diff[j]).real();
            out[i * n_eps2 + j] = std::isfinite(val) ? val : 0.0;
        }
    }
}

void regularized_J(double E1, const double* E2, int n_E2,
                   const double* epsilon, int n_eps,
                   double T, double* out) {
    if (n_E2 <= 0 || n_eps <= 0) {
        return;
    }

    double beta = 1.0 / (KB_EV * T);
    double coeff = beta / (2.0 * PI);

    std::vector<double> bose_vals(n_E2);
    std::vector<Complex> trig_a2(n_eps);

    for (int j = 0; j < n_E2; ++j) {
        bose_vals[j] = bose_fcn(E2[j] - E1, T);
    }

    for (int i = 0; i < n_eps; ++i) {
        Complex a2 = 0.5 + Complex(0.0, 1.0) * coeff * (E1 - epsilon[i]);
        trig_a2[i] = trigamma(a2);
    }

    for (int i = 0; i < n_eps; ++i) {
        for (int j = 0; j < n_E2; ++j) {
            Complex a1 = 0.5 + Complex(0.0, 1.0) * coeff * (E2[j] - epsilon[i]);
            double val = coeff * bose_vals[j] * (trigamma(a1) - trig_a2[i]).imag();
            out[i * n_E2 + j] = std::isfinite(val) ? val : 0.0;
        }
    }
}

void regularized_J_matrix(double E1, const double* E2, int n_E2,
                          const double* epsilon_matrix, int eps_rows, int eps_cols,
                          double T, double* out) {
    if (n_E2 <= 0 || eps_rows <= 0 || eps_cols <= 0) {
        return;
    }

    double beta = 1.0 / (KB_EV * T);
    double coeff = beta / (2.0 * PI);

    std::vector<double> bose_vals(n_E2);
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
            Complex a1 = 0.5 + Complex(0.0, 1.0) * coeff * (E2[j] - eps);
            Complex a2 = 0.5 + Complex(0.0, 1.0) * coeff * (E1 - eps);

            double val = coeff * bose_vals[j] * (trigamma(a1) - trigamma(a2)).imag();
            out[i * n_E2 + j] = std::isfinite(val) ? val : 0.0;
        }
    }
}

}
