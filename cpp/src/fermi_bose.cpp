#include "fermi_bose.hpp"

#include <cmath>
#include "constants.hpp"

namespace fc {

double fermi(double x, double mu, double T) {
    return 1.0 / (std::exp((x - mu) / (KB_EV * T)) + 1.0);
}

double bose_fcn(double x, double T) {
    double beta = 1.0 / (KB_EV * T);
    return 1.0 / (std::exp(x * beta) - 1.0);
}

}
