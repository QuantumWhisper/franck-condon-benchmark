#include "fc_matrix.h"

#include <math.h>

#include "laguerre.h"

void fc_cache_init(FCCache *fc, double lambda) {
    if (fc == 0) {
        return;
    }
    fc->lambda = lambda;
    for (int i = 0; i < FC_MAX_N; ++i) {
        for (int j = 0; j < FC_MAX_N; ++j) {
            fc->valid[i][j] = 0;
            fc->cache[i][j] = 0.0;
        }
    }
}

double fc_matrix_single(int q1, int q2, double lambda) {
    int q = (q1 < q2) ? q1 : q2;
    int Q = (q1 > q2) ? q1 : q2;

    if (q1 == q2 && q1 == 0) {
        return exp(-lambda * lambda / 2.0);
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
        (double)(Q - q) * log(lambda) -
        lambda * lambda / 2.0 +
        0.5 * (lgamma((double)q + 1.0) - lgamma((double)Q + 1.0));

    double M = sign_factor * exp(log_coeff) * L;
    if (isnan(M) || isinf(M)) {
        return 0.0;
    }
    return M;
}

double fc_cache_get(FCCache *fc, int q1, int q2) {
    if (fc == 0 || q1 < 0 || q2 < 0 || q1 >= FC_MAX_N || q2 >= FC_MAX_N) {
        return 0.0;
    }

    if (!fc->valid[q1][q2]) {
        fc->cache[q1][q2] = fc_matrix_single(q1, q2, fc->lambda);
        fc->valid[q1][q2] = 1;
    }
    return fc->cache[q1][q2];
}

void fc_matrix_row(FCCache *fc, int q1, const int *q2_range, int nq2, double *out) {
    if (fc == 0 || q2_range == 0 || out == 0 || nq2 <= 0) {
        return;
    }

    for (int i = 0; i < nq2; ++i) {
        out[i] = fc_cache_get(fc, q1, q2_range[i]);
    }
}
