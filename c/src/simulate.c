#include "simulate.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "constants.h"
#include "current.h"
#include "fc_matrix.h"
#include "rate.h"

void simulation_result_free(SimulationResult *res) {
    if (res == NULL) {
        return;
    }
    free(res->Vsd);
    free(res->I_tol);
    free(res->I_seq);
    free(res->I_cot);
    res->Vsd = NULL;
    res->I_tol = NULL;
    res->I_seq = NULL;
    res->I_cot = NULL;
    res->nVsd = 0;
}

SimulationResult simulate_iv(int N, double vmode, double alphaL, double alphaR,
                             double lambda, const double *Vsd_vec, int nVsd,
                             double T, double eta, double Vg, double tau,
                             int verbose) {
    SimulationResult res;
    res.Vsd = NULL;
    res.I_tol = NULL;
    res.I_seq = NULL;
    res.I_cot = NULL;
    res.nVsd = nVsd;

    if (Vsd_vec == NULL || nVsd <= 0 || N <= 0) {
        return res;
    }

    res.Vsd = (double *)malloc((size_t)nVsd * sizeof(double));
    res.I_tol = (double *)malloc((size_t)nVsd * sizeof(double));
    res.I_seq = (double *)malloc((size_t)nVsd * sizeof(double));
    res.I_cot = (double *)malloc((size_t)nVsd * sizeof(double));
    if (res.Vsd == NULL || res.I_tol == NULL || res.I_seq == NULL || res.I_cot == NULL) {
        simulation_result_free(&res);
        return res;
    }
    memcpy(res.Vsd, Vsd_vec, (size_t)nVsd * sizeof(double));

    FCCache *fc = (FCCache *)malloc(sizeof(FCCache));
    if (fc == NULL) {
        simulation_result_free(&res);
        return res;
    }
    fc_cache_init(fc, lambda);
    fc_cache_populate(fc, FC_MAX_N);

    #pragma omp parallel for schedule(dynamic, 4)
    for (int vv = 0; vv < nVsd; ++vv) {
        const double v = Vsd_vec[vv];
        if (verbose) {
#ifdef _OPENMP
            if (omp_get_thread_num() == 0)
#endif
            {
                fprintf(stderr, "\rBias point %d/%d (Vsd = %.4f V)", vv + 1, nVsd, v);
                fflush(stderr);
            }
        }

        RateStore store;
        rate_store_init(&store, N);

        calculate_all_rateW(&store, N, vmode, alphaL, alphaR, lambda, v, T, eta, 1, Vg, fc);
        calculate_all_rateW(&store, N, vmode, alphaL, alphaR, lambda, v, T, eta, -1, Vg, fc);

        CurrentResult cr =
            current_from_rate_equations(v, N, vmode, alphaL, alphaR, lambda, T, eta, Vg, tau,
                                        &store);

        if (v != 0.0) {
            const double s = (v > 0.0) ? -1.0 : 1.0;
            cr.I_tol *= s;
            cr.I_seq *= s;
            cr.I_cot *= s;
        }

        res.I_tol[vv] = cr.I_tol * ELEMENTARY_CHARGE;
        res.I_seq[vv] = cr.I_seq * ELEMENTARY_CHARGE;
        res.I_cot[vv] = cr.I_cot * ELEMENTARY_CHARGE;

        rate_store_free(&store);
    }

    if (verbose) {
        fprintf(stderr, "\n");
    }

    free(fc);
    return res;
}
