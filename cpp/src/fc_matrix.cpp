#include "fc_matrix.hpp"

#include <cmath>
#include "laguerre.hpp"

namespace fc {

FCCache::FCCache(double lambda) : lambda(lambda) {
    for (int i = 0; i < FC_MAX_N; ++i) {
        for (int j = 0; j < FC_MAX_N; ++j) {
            valid[i][j] = false;
            cache[i][j] = 0.0;
        }
    }
}

double fc_matrix_single(int q1, int q2, double lambda) {
    int q = (q1 < q2) ? q1 : q2;
    int Q = (q1 > q2) ? q1 : q2;

    if (q1 == q2 && q1 == 0) {
        return std::exp(-lambda * lambda / 2.0);
    }

    double sign_factor = 1.0;
    if (q2 < q1 && ((q1 - q2) % 2 != 0)) {
        sign_factor = -1.0;
    }

    double L = laguerre_L(q, Q - q, lambda * lambda);

    if (lambda == 0.0) {
        if (q1 == q2) {
            return L;
        }
        return 0.0;
    }

    double log_coeff =
        static_cast<double>(Q - q) * std::log(lambda) -
        lambda * lambda / 2.0 +
        0.5 * (std::lgamma(static_cast<double>(q) + 1.0) -
               std::lgamma(static_cast<double>(Q) + 1.0));

    double M = sign_factor * std::exp(log_coeff) * L;
    if (std::isnan(M) || std::isinf(M)) {
        return 0.0;
    }
    return M;
}

double fc_cache_get(FCCache& fc, int q1, int q2) {
    if (q1 < 0 || q2 < 0 || q1 >= FC_MAX_N || q2 >= FC_MAX_N) {
        return 0.0;
    }

    if (!fc.valid[q1][q2]) {
        fc.cache[q1][q2] = fc_matrix_single(q1, q2, fc.lambda);
        fc.valid[q1][q2] = true;
    }
    return fc.cache[q1][q2];
}

void fc_matrix_row(FCCache& fc, int q1, const int* q2_range, int nq2, double* out) {
    if (nq2 <= 0) {
        return;
    }

    for (int i = 0; i < nq2; ++i) {
        out[i] = fc_cache_get(fc, q1, q2_range[i]);
    }
}

void fc_cache_populate(FCCache& fc, int N) {
    for (int q1 = 0; q1 < N; ++q1) {
        for (int q2 = 0; q2 < N; ++q2) {
            fc_cache_get(fc, q1, q2);
        }
    }
}

}
