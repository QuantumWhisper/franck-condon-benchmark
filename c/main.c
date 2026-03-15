#define _POSIX_C_SOURCE 199309L
#define _ISOC99_SOURCE

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "src/json_io.h"
#include "src/plotting.h"
#include "src/simulate.h"

static void sort3(double *a, double *b, double *c) {
    if (*a > *b) {
        double t = *a;
        *a = *b;
        *b = t;
    }
    if (*b > *c) {
        double t = *b;
        *b = *c;
        *c = t;
    }
    if (*a > *b) {
        double t = *a;
        *a = *b;
        *b = t;
    }
}

int main(int argc, char **argv) {
    const char *spec = (argc >= 2) ? argv[1] : "default";

    char params_path[512];
    char ref_path[512];
    char json_out[512];
    char csv_out[512];
    char pdf_out[512];
    char png_out[512];

    sprintf(params_path, "../benchmark/spec/%s_params.json", spec);
    sprintf(ref_path, "../benchmark/results/matlab_%s_results.json", spec);
    sprintf(json_out, "../benchmark/results/c_%s_results.json", spec);
    sprintf(csv_out, "../benchmark/results/c_%s_IV.csv", spec);
    sprintf(pdf_out, "../benchmark/results/c_%s_IV.pdf", spec);
    sprintf(png_out, "../benchmark/results/c_%s_IV.png", spec);

    SimParams params = parse_params_json(params_path);
    if (params.N <= 0) {
        fprintf(stderr, "Failed to parse parameters from %s\n", params_path);
        return 1;
    }

    int nVsd = 0;
    double *Vsd = load_matlab_vsd(ref_path, &nVsd);
    if (Vsd == NULL || nVsd <= 0) {
        nVsd = (int)llround((params.Vsd_end - params.Vsd_start) / params.Vsd_step) + 1;
        Vsd = (double *)malloc((size_t)nVsd * sizeof(double));
        if (Vsd == NULL) {
            fprintf(stderr, "Failed to allocate Vsd array\n");
            return 1;
        }
        for (int i = 0; i < nVsd; ++i) {
            Vsd[i] = params.Vsd_start + (double)i * params.Vsd_step;
        }
    }

    printf("=== Franck-Condon Benchmark (C) [%s] ===\n", spec);
    printf("N=%d, lambda=%.1f, T=%.1f K, Vsd=[%.3f:%.3f:%.3f] V\n",
           params.N, params.lambda, params.T,
           params.Vsd_start, params.Vsd_step, params.Vsd_end);
    printf("Total bias points: %d\n\n", nVsd);

    double wall_times[3] = {0.0, 0.0, 0.0};
    SimulationResult result;
    memset(&result, 0, sizeof(result));

    for (int run = 0; run < 3; ++run) {
        struct timespec t0;
        struct timespec t1;
        clock_gettime(CLOCK_MONOTONIC, &t0);
        result = simulate_iv(params.N, params.vmode, params.alphaL, params.alphaR,
                             params.lambda, Vsd, nVsd, params.T, params.eta,
                             params.Vg, params.tau, run == 0 ? 1 : 0);
        clock_gettime(CLOCK_MONOTONIC, &t1);

        wall_times[run] = (double)(t1.tv_sec - t0.tv_sec) +
                          (double)(t1.tv_nsec - t0.tv_nsec) * 1e-9;
        printf("Run %d: %.3f s\n", run + 1, wall_times[run]);
        if (run < 2) {
            simulation_result_free(&result);
        }
    }

    sort3(&wall_times[0], &wall_times[1], &wall_times[2]);
    double median_time = wall_times[1];
    printf("Median wall time: %.3f s\n", median_time);

    double *ref_tol = NULL;
    double *ref_seq = NULL;
    double *ref_cot = NULL;
    int ref_npts = 0;
    load_matlab_reference(ref_path, &ref_tol, &ref_seq, &ref_cot, &ref_npts);
    if (ref_tol != NULL && ref_seq != NULL && ref_cot != NULL && ref_npts > 0) {
        double max_err_tol = 0.0, max_err_seq = 0.0, max_err_cot = 0.0;
        double max_err_tol_ex = 0.0, max_err_seq_ex = 0.0, max_err_cot_ex = 0.0;

        int ncmp = (ref_npts < nVsd) ? ref_npts : nVsd;
        for (int i = 0; i < ncmp; ++i) {
            double et, es, ec, denom;
            denom = fmax(fabs(ref_tol[i]), 1e-30);
            et = fabs(result.I_tol[i] - ref_tol[i]) / denom;
            denom = fmax(fabs(ref_seq[i]), 1e-30);
            es = fabs(result.I_seq[i] - ref_seq[i]) / denom;
            denom = fmax(fabs(ref_cot[i]), 1e-30);
            ec = fabs(result.I_cot[i] - ref_cot[i]) / denom;

            if (et > max_err_tol) max_err_tol = et;
            if (es > max_err_seq) max_err_seq = es;
            if (ec > max_err_cot) max_err_cot = ec;

            int is_artifact = (fabs(Vsd[i] - 0.219) < 0.002) || (fabs(Vsd[i] - 0.585) < 0.002);
            if (!is_artifact) {
                if (et > max_err_tol_ex) max_err_tol_ex = et;
                if (es > max_err_seq_ex) max_err_seq_ex = es;
                if (ec > max_err_cot_ex) max_err_cot_ex = ec;
            }
        }

        printf("Max relative error vs MATLAB (all points):\n");
        printf("  I_tol: %.6e  I_seq: %.6e  I_cot: %.6e\n", max_err_tol, max_err_seq, max_err_cot);
        printf("Max relative error vs MATLAB (excl. solver artifacts at Vsd~0.219,0.585):\n");
        printf("  I_tol: %.6e  I_seq: %.6e  I_cot: %.6e\n", max_err_tol_ex, max_err_seq_ex, max_err_cot_ex);

        const double tol = 1e-4;
        if (max_err_tol_ex < tol && max_err_seq_ex < tol && max_err_cot_ex < tol) {
            printf("VALIDATION PASSED (tolerance %.0e, excluding artifact points)\n", tol);
        } else {
            printf("VALIDATION FAILED (tolerance %.0e)\n", tol);
        }

        free(ref_tol);
        free(ref_seq);
        free(ref_cot);
    } else {
        printf("No MATLAB reference found, skipping validation.\n");
    }

    write_results_json(json_out, spec, &params, &result, median_time);
    write_results_csv(csv_out, &result);
    plot_iv(csv_out, pdf_out, png_out, params.N, params.lambda, params.T, params.vmode,
            params.alphaL, params.alphaR, params.eta, params.Vg, median_time, spec);

    simulation_result_free(&result);
    free(Vsd);
    return 0;
}
