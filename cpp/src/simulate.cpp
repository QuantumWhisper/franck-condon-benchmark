#include "simulate.hpp"

#include <cstdio>
#include <cstring>
#include <omp.h>

#include "constants.hpp"
#include "current.hpp"
#include "fc_matrix.hpp"
#include "rate.hpp"

namespace fc {

SimulationResult simulate_iv(int N, double vmode, double alphaL, double alphaR,
                             double lambda, const double* Vsd_vec, int nVsd,
                             double T, double eta, double Vg, double tau,
                             int verbose) {
    SimulationResult res;

    if (Vsd_vec == nullptr || nVsd <= 0 || N <= 0) {
        return res;
    }

    res.Vsd.assign(Vsd_vec, Vsd_vec + nVsd);
    res.I_tol.resize(nVsd, 0.0);
    res.I_seq.resize(nVsd, 0.0);
    res.I_cot.resize(nVsd, 0.0);

    FCCache fc(lambda);
    fc_cache_populate(fc, N);

    #pragma omp parallel for schedule(dynamic, 4)
    for (int vv = 0; vv < nVsd; ++vv) {
        const double v = Vsd_vec[vv];

        if (verbose && omp_get_thread_num() == 0) {
            std::fprintf(stderr, "\rBias point %d/%d (Vsd = %.4f V)", vv + 1, nVsd, v);
            std::fflush(stderr);
        }

        RateStore store(N);

        FCCache fc_local(lambda);
        fc_cache_populate(fc_local, N);

        calculate_all_rateW(store, N, vmode, alphaL, alphaR, lambda, v, T, eta, 1, Vg, fc_local);
        calculate_all_rateW(store, N, vmode, alphaL, alphaR, lambda, v, T, eta, -1, Vg, fc_local);

        CurrentResult cr =
            current_from_rate_equations(v, N, vmode, alphaL, alphaR, lambda, T, eta, Vg, tau,
                                        store);

        if (v != 0.0) {
            const double s = (v > 0.0) ? -1.0 : 1.0;
            cr.I_tol *= s;
            cr.I_seq *= s;
            cr.I_cot *= s;
        }

        res.I_tol[vv] = cr.I_tol * ELEMENTARY_CHARGE;
        res.I_seq[vv] = cr.I_seq * ELEMENTARY_CHARGE;
        res.I_cot[vv] = cr.I_cot * ELEMENTARY_CHARGE;
    }

    if (verbose) {
        std::fprintf(stderr, "\n");
    }

    return res;
}

}
