#pragma once

#include <string>
#include <vector>
#include "simulate.hpp"

namespace fc {

struct SimParams {
    int N = 0;
    double vmode = 0.0;
    double alphaL = 0.0;
    double alphaR = 0.0;
    double lambda = 0.0;
    double T = 0.0;
    double eta = 0.0;
    double Vg = 0.0;
    double tau = 0.0;
    double Vsd_start = 0.0;
    double Vsd_end = 0.0;
    double Vsd_step = 0.0;
};

struct MatlabReference {
    std::vector<double> I_tol;
    std::vector<double> I_seq;
    std::vector<double> I_cot;
};

SimParams parse_params_json(const std::string& filepath);
std::vector<double> load_matlab_vsd(const std::string& filepath);
MatlabReference load_matlab_reference(const std::string& filepath);
void write_results_json(const std::string& filepath, const std::string& spec,
                        const SimParams& params, const SimulationResult& res,
                        double wall_time);
void write_results_csv(const std::string& filepath, const SimulationResult& res);

}
