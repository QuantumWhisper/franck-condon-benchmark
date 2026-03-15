#ifndef FC_SIMULATE_H
#define FC_SIMULATE_H

typedef struct {
    double *Vsd;
    double *I_tol;
    double *I_seq;
    double *I_cot;
    int nVsd;
} SimulationResult;

SimulationResult simulate_iv(int N, double vmode, double alphaL, double alphaR,
                             double lambda, const double *Vsd_vec, int nVsd,
                             double T, double eta, double Vg, double tau,
                             int verbose);
void simulation_result_free(SimulationResult *res);

#endif
