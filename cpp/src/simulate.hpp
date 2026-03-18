#pragma once

#include <vector>

namespace fc {

struct SimulationResult {
    std::vector<double> Vsd;
    std::vector<double> I_tol;
    std::vector<double> I_seq;
    std::vector<double> I_cot;
};

SimulationResult simulate_iv(int N, double vmode, double alphaL, double alphaR,
                             double lambda, const double* Vsd_vec, int nVsd,
                             double T, double eta, double Vg, double tau,
                             int verbose);

}
