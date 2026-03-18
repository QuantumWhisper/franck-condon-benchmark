#pragma once

#include "rate.hpp"

namespace fc {

struct CurrentResult {
    double I_tol = 0.0;
    double I_seq = 0.0;
    double I_cot = 0.0;
};

CurrentResult current_from_rate_equations(
    double Vsd, int N, double vmode, double alphaL, double alphaR,
    double lambda, double T, double eta, double Vg, double tau,
    const RateStore& store);

}
