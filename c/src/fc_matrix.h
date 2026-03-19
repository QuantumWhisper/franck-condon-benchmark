#ifndef FC_FC_MATRIX_H
#define FC_FC_MATRIX_H

#define FC_MAX_N 256

typedef struct {
    double cache[FC_MAX_N][FC_MAX_N];
    int valid[FC_MAX_N][FC_MAX_N];
    double lambda;
} FCCache;

void fc_cache_init(FCCache *fc, double lambda);
double fc_matrix_single(int q1, int q2, double lambda);
double fc_cache_get(FCCache *fc, int q1, int q2);
void fc_matrix_row(FCCache *fc, int q1, const int *q2_range, int nq2, double *out);
void fc_cache_populate(FCCache *fc, int max_q);

#endif
