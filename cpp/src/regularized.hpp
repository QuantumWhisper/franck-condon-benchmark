#pragma once

namespace fc {

void regularized_I(double E1, double E2,
                   const double* epsilon1, int n_eps1,
                   const double* epsilon2, int n_eps2,
                   double T, double* out);

void regularized_J(double E1, const double* E2, int n_E2,
                   const double* epsilon, int n_eps,
                   double T, double* out);

void regularized_J_matrix(double E1, const double* E2, int n_E2,
                          const double* epsilon_matrix, int eps_rows, int eps_cols,
                          double T, double* out);

}
