#pragma once

namespace fc {

inline constexpr int FC_MAX_N = 256;

struct FCCache {
    double cache[FC_MAX_N][FC_MAX_N];
    bool valid[FC_MAX_N][FC_MAX_N];
    double lambda;

    explicit FCCache(double lambda);
};

double fc_matrix_single(int q1, int q2, double lambda);
double fc_cache_get(FCCache& fc, int q1, int q2);
void fc_matrix_row(FCCache& fc, int q1, const int* q2_range, int nq2, double* out);
void fc_cache_populate(FCCache& fc, int N);

}
