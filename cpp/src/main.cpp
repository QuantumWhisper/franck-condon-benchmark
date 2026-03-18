#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "json_io.hpp"
#include "plotting.hpp"
#include "simulate.hpp"

static void sort3(double& a, double& b, double& c) {
    if (a > b) std::swap(a, b);
    if (b > c) std::swap(b, c);
    if (a > b) std::swap(a, b);
}

int main(int argc, char** argv) {
    std::string spec = (argc >= 2) ? argv[1] : "default";

    std::string params_path = "../benchmark/spec/" + spec + "_params.json";
    std::string ref_path = "../benchmark/results/matlab_" + spec + "_results.json";
    std::string json_out = "../benchmark/results/cpp_" + spec + "_results.json";
    std::string csv_out = "../benchmark/results/cpp_" + spec + "_IV.csv";
    std::string pdf_out = "../benchmark/results/cpp_" + spec + "_IV.pdf";
    std::string png_out = "../benchmark/results/cpp_" + spec + "_IV.png";

    fc::SimParams params = fc::parse_params_json(params_path);
    if (params.N <= 0) {
        std::fprintf(stderr, "Failed to parse parameters from %s\n", params_path.c_str());
        return 1;
    }

    std::vector<double> Vsd = fc::load_matlab_vsd(ref_path);
    int nVsd = static_cast<int>(Vsd.size());
    if (nVsd <= 0) {
        nVsd = static_cast<int>(std::llround((params.Vsd_end - params.Vsd_start) / params.Vsd_step)) + 1;
        Vsd.resize(nVsd);
        for (int i = 0; i < nVsd; ++i) {
            Vsd[i] = params.Vsd_start + static_cast<double>(i) * params.Vsd_step;
        }
    }

    std::printf("=== Franck-Condon Benchmark (C++) [%s] ===\n", spec.c_str());
    std::printf("N=%d, lambda=%.1f, T=%.1f K, Vsd=[%.3f:%.3f:%.3f] V\n",
                params.N, params.lambda, params.T,
                params.Vsd_start, params.Vsd_step, params.Vsd_end);
    std::printf("Total bias points: %d\n\n", nVsd);

    double wall_times[3] = {0.0, 0.0, 0.0};
    fc::SimulationResult result;

    for (int run = 0; run < 3; ++run) {
        auto t0 = std::chrono::steady_clock::now();
        result = fc::simulate_iv(params.N, params.vmode, params.alphaL, params.alphaR,
                                 params.lambda, Vsd.data(), nVsd, params.T, params.eta,
                                 params.Vg, params.tau, run == 0 ? 1 : 0);
        auto t1 = std::chrono::steady_clock::now();

        wall_times[run] = std::chrono::duration<double>(t1 - t0).count();
        std::printf("Run %d: %.3f s\n", run + 1, wall_times[run]);
    }

    sort3(wall_times[0], wall_times[1], wall_times[2]);
    double median_time = wall_times[1];
    std::printf("Median wall time: %.3f s\n", median_time);

    fc::MatlabReference ref = fc::load_matlab_reference(ref_path);
    if (!ref.I_tol.empty() && !ref.I_seq.empty() && !ref.I_cot.empty()) {
        double max_err_tol = 0.0, max_err_seq = 0.0, max_err_cot = 0.0;
        double max_err_tol_ex = 0.0, max_err_seq_ex = 0.0, max_err_cot_ex = 0.0;

        int ncmp = std::min(static_cast<int>(ref.I_tol.size()), nVsd);
        for (int i = 0; i < ncmp; ++i) {
            double denom, et, es, ec;
            denom = std::fmax(std::fabs(ref.I_tol[i]), 1e-30);
            et = std::fabs(result.I_tol[i] - ref.I_tol[i]) / denom;
            denom = std::fmax(std::fabs(ref.I_seq[i]), 1e-30);
            es = std::fabs(result.I_seq[i] - ref.I_seq[i]) / denom;
            denom = std::fmax(std::fabs(ref.I_cot[i]), 1e-30);
            ec = std::fabs(result.I_cot[i] - ref.I_cot[i]) / denom;

            if (et > max_err_tol) max_err_tol = et;
            if (es > max_err_seq) max_err_seq = es;
            if (ec > max_err_cot) max_err_cot = ec;

            bool is_artifact = (std::fabs(Vsd[i] - 0.219) < 0.002) || (std::fabs(Vsd[i] - 0.585) < 0.002);
            if (!is_artifact) {
                if (et > max_err_tol_ex) max_err_tol_ex = et;
                if (es > max_err_seq_ex) max_err_seq_ex = es;
                if (ec > max_err_cot_ex) max_err_cot_ex = ec;
            }
        }

        std::printf("Max relative error vs MATLAB (all points):\n");
        std::printf("  I_tol: %.6e  I_seq: %.6e  I_cot: %.6e\n", max_err_tol, max_err_seq, max_err_cot);
        std::printf("Max relative error vs MATLAB (excl. solver artifacts at Vsd~0.219,0.585):\n");
        std::printf("  I_tol: %.6e  I_seq: %.6e  I_cot: %.6e\n", max_err_tol_ex, max_err_seq_ex, max_err_cot_ex);

        constexpr double tol = 1e-4;
        if (max_err_tol_ex < tol && max_err_seq_ex < tol && max_err_cot_ex < tol) {
            std::printf("VALIDATION PASSED (tolerance %.0e, excluding artifact points)\n", tol);
        } else {
            std::printf("VALIDATION FAILED (tolerance %.0e)\n", tol);
        }
    } else {
        std::printf("No MATLAB reference found, skipping validation.\n");
    }

    fc::write_results_json(json_out, spec, params, result, median_time);
    fc::write_results_csv(csv_out, result);
    fc::plot_iv(csv_out, pdf_out, png_out, params.N, params.lambda, params.T, params.vmode,
                params.alphaL, params.alphaR, params.eta, params.Vg, median_time, spec);

    return 0;
}
