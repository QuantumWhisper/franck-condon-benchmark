#include "solver.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>

typedef struct gsl_matrix gsl_matrix;
typedef struct gsl_vector gsl_vector;

gsl_matrix *gsl_matrix_alloc(size_t n1, size_t n2);
void gsl_matrix_free(gsl_matrix *m);
void gsl_matrix_set(gsl_matrix *m, size_t i, size_t j, double x);

gsl_vector *gsl_vector_alloc(size_t n);
gsl_vector *gsl_vector_calloc(size_t n);
void gsl_vector_free(gsl_vector *v);
void gsl_vector_set(gsl_vector *v, size_t i, double x);
double gsl_vector_get(const gsl_vector *v, size_t i);

int gsl_linalg_QR_decomp(gsl_matrix *A, gsl_vector *tau);
int gsl_linalg_QR_solve(const gsl_matrix *QR,
                        const gsl_vector *tau,
                        const gsl_vector *b,
                        gsl_vector *x);

void solve_steady_state(double *W, double *P, int N) {
    int dim = 2 * N;
    if (W == NULL || P == NULL || dim <= 0) {
        return;
    }

    gsl_matrix *A = gsl_matrix_alloc((size_t)dim, (size_t)dim);
    gsl_vector *d = gsl_vector_calloc((size_t)dim);
    gsl_vector *tau = gsl_vector_alloc((size_t)dim);
    gsl_vector *x = gsl_vector_alloc((size_t)dim);

    if (A == NULL || d == NULL || tau == NULL || x == NULL) {
        gsl_matrix_free(A);
        gsl_vector_free(d);
        gsl_vector_free(tau);
        gsl_vector_free(x);
        return;
    }

    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            double val = W[i * dim + j];
            if (i == dim - 1) {
                val = 1.0;
            }
            gsl_matrix_set(A, (size_t)i, (size_t)j, val);
        }
    }
    gsl_vector_set(d, (size_t)(dim - 1), 1.0);

    gsl_linalg_QR_decomp(A, tau);
    gsl_linalg_QR_solve(A, tau, d, x);

    double sum = 0.0;
    for (int i = 0; i < dim; ++i) {
        double v = gsl_vector_get(x, (size_t)i);
        if (!isfinite(v) || v < 0.0) {
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

    gsl_matrix_free(A);
    gsl_vector_free(d);
    gsl_vector_free(tau);
    gsl_vector_free(x);
}
