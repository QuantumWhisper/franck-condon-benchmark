#ifndef FC_JSON_IO_H
#define FC_JSON_IO_H

#include "simulate.h"

#include "../cJSON.h"

typedef struct {
    int N;
    double vmode;
    double alphaL;
    double alphaR;
    double lambda;
    double T;
    double eta;
    double Vg;
    double tau;
    double Vsd_start;
    double Vsd_end;
    double Vsd_step;
} SimParams;

SimParams parse_params_json(const char *filepath);
double *load_matlab_vsd(const char *filepath, int *nVsd);
void load_matlab_reference(const char *filepath,
                           double **ref_I_tol, double **ref_I_seq,
                           double **ref_I_cot, int *npts);
void write_results_json(const char *filepath, const char *spec,
                        const SimParams *params, const SimulationResult *res,
                        double wall_time);
void write_results_csv(const char *filepath, const SimulationResult *res);

#endif
