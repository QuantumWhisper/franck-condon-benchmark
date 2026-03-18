#include "solver.hpp"

#include <cmath>
#include <Eigen/Dense>

namespace fc {

void solve_steady_state(double* W, double* P, int N) {
    int dim = 2 * N;
    if (dim <= 0) {
        return;
    }

    Eigen::MatrixXd A(dim, dim);
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            A(i, j) = W[i * dim + j];
        }
    }

    A.row(dim - 1).setOnes();

    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(dim);
    rhs(dim - 1) = 1.0;

    Eigen::VectorXd x = A.householderQr().solve(rhs);

    double sum = 0.0;
    for (int i = 0; i < dim; ++i) {
        double v = x(i);
        if (!std::isfinite(v) || v < 0.0) {
            v = 0.0;
        }
        P[i] = v;
        sum += v;
    }

    if (sum > 0.0) {
        for (int i = 0; i < dim; ++i) {
            P[i] /= sum;
        }
    } else {
        for (int i = 0; i < dim; ++i) {
            P[i] = 0.0;
        }
        P[dim - 1] = 1.0;
    }
}

}
