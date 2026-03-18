#include "current.hpp"

#include <vector>

#include "matrix.hpp"
#include "solver.hpp"

namespace fc {

CurrentResult current_from_rate_equations(
    double Vsd, int N, double vmode, double alphaL, double alphaR,
    double lambda, double T, double eta, double Vg, double tau,
    const RateStore& store) {
    (void)Vsd;
    (void)alphaL;
    (void)alphaR;
    (void)lambda;
    (void)eta;
    (void)Vg;

    CurrentResult out;

    if (N <= 0) {
        return out;
    }

    const int dim = 2 * N;
    std::vector<double> W(dim * dim);
    std::vector<double> P(dim);

    generate_matrix_W(store, N, vmode, T, tau, W.data());
    solve_steady_state(W.data(), P.data(), N);

    constexpr int leadR = -1;
    constexpr int leadL = 1;

    double I_seq0 = 0.0;
    for (int q1 = 0; q1 < N; ++q1) {
        double sum_diff = 0.0;
        for (int q2 = 0; q2 < N; ++q2) {
            const double w_R = rateW_from_store(store, 0, 1, q1, q2, leadR);
            const double w_L = rateW_from_store(store, 0, 1, q1, q2, leadL);
            sum_diff += (w_R - w_L);
        }
        I_seq0 += P[q1] * sum_diff;
    }

    double I_seq1 = 0.0;
    for (int q1 = 0; q1 < N; ++q1) {
        double sum_diff = 0.0;
        for (int q2 = 0; q2 < N; ++q2) {
            const double w_R = rateW_from_store(store, 1, 0, q1, q2, leadR);
            const double w_L = rateW_from_store(store, 1, 0, q1, q2, leadL);
            sum_diff += (w_L - w_R);
        }
        I_seq1 += P[N + q1] * sum_diff;
    }

    const double I_seq = I_seq0 + I_seq1;

    double I_cot = 0.0;
    for (int nn = 0; nn <= 1; ++nn) {
        for (int q1 = 0; q1 < N; ++q1) {
            double sum_diff = 0.0;
            for (int q2 = 0; q2 < N; ++q2) {
                const double w_RL = rateW_from_store(store, nn, nn, q1, q2, leadR);
                const double w_LR = rateW_from_store(store, nn, nn, q1, q2, leadL);
                sum_diff += (w_RL - w_LR);
            }
            I_cot += P[nn * N + q1] * sum_diff;
        }
    }

    out.I_seq = I_seq;
    out.I_cot = I_cot;
    out.I_tol = I_seq + I_cot;

    return out;
}

}
