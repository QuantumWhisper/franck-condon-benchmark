#include "json_io.hpp"

#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <limits>

#include <nlohmann/json.hpp>

namespace fc {

SimParams parse_params_json(const std::string& filepath) {
    SimParams p;
    p.tau = std::numeric_limits<double>::infinity();

    std::ifstream f(filepath);
    if (!f.is_open()) {
        return p;
    }

    nlohmann::json root;
    try {
        root = nlohmann::json::parse(f);
    } catch (...) {
        return p;
    }

    auto& params = root["parameters"];
    auto& sweep = root["bias_sweep"];

    p.N = params.value("N", 0);
    p.vmode = params.value("vmode", 0.0);
    p.alphaL = params.value("alphaL", 0.0);
    p.alphaR = params.value("alphaR", 0.0);
    p.lambda = params.value("lambda", 0.0);
    p.T = params.value("T", 0.0);
    p.eta = params.value("eta", 0.0);
    p.Vg = params.value("Vg", 0.0);

    auto tau_val = params["tau"];
    if (tau_val.is_string() && tau_val.get<std::string>() == "Inf") {
        p.tau = std::numeric_limits<double>::infinity();
    } else if (tau_val.is_number()) {
        p.tau = tau_val.get<double>();
    }

    p.Vsd_start = sweep.value("Vsd_start", 0.0);
    p.Vsd_end = sweep.value("Vsd_end", 0.0);
    p.Vsd_step = sweep.value("Vsd_step", 0.0);

    return p;
}

std::vector<double> load_matlab_vsd(const std::string& filepath) {
    std::ifstream f(filepath);
    if (!f.is_open()) {
        return {};
    }

    nlohmann::json root;
    try {
        root = nlohmann::json::parse(f);
    } catch (...) {
        return {};
    }

    if (!root.contains("Vsd") || !root["Vsd"].is_array()) {
        return {};
    }

    return root["Vsd"].get<std::vector<double>>();
}

MatlabReference load_matlab_reference(const std::string& filepath) {
    MatlabReference ref;

    std::ifstream f(filepath);
    if (!f.is_open()) {
        return ref;
    }

    nlohmann::json root;
    try {
        root = nlohmann::json::parse(f);
    } catch (...) {
        return ref;
    }

    if (root.contains("I_tol") && root["I_tol"].is_array()) {
        ref.I_tol = root["I_tol"].get<std::vector<double>>();
    }
    if (root.contains("I_seq") && root["I_seq"].is_array()) {
        ref.I_seq = root["I_seq"].get<std::vector<double>>();
    }
    if (root.contains("I_cot") && root["I_cot"].is_array()) {
        ref.I_cot = root["I_cot"].get<std::vector<double>>();
    }

    return ref;
}

void write_results_json(const std::string& filepath, const std::string& spec,
                        const SimParams& params, const SimulationResult& res,
                        double wall_time) {
    FILE* fp = std::fopen(filepath.c_str(), "w");
    if (fp == nullptr) {
        return;
    }

    std::time_t now = std::time(nullptr);
    std::tm* lt = std::localtime(&now);
    char tbuf[64] = {0};
    if (lt != nullptr) {
        std::strftime(tbuf, sizeof(tbuf), "%Y-%m-%d %H:%M:%S", lt);
    }

    int n = static_cast<int>(res.Vsd.size());

    std::fprintf(fp, "{\n");
    std::fprintf(fp, "  \"spec\": \"%s\",\n", spec.c_str());
    std::fprintf(fp, "  \"language\": \"C++\",\n");
#if defined(__clang__)
    std::fprintf(fp, "  \"version\": \"clang %s\",\n", __clang_version__);
#elif defined(__GNUC__)
    std::fprintf(fp, "  \"version\": \"gcc %s\",\n", __VERSION__);
#else
    std::fprintf(fp, "  \"version\": \"unknown\",\n");
#endif
    std::fprintf(fp, "  \"wall_time_seconds\": %.17g,\n", wall_time);
    std::fprintf(fp, "  \"parameters\": {\n");
    std::fprintf(fp, "    \"N\": %d,\n", params.N);
    std::fprintf(fp, "    \"vmode\": %.17g,\n", params.vmode);
    std::fprintf(fp, "    \"alphaL\": %.17g,\n", params.alphaL);
    std::fprintf(fp, "    \"alphaR\": %.17g,\n", params.alphaR);
    std::fprintf(fp, "    \"lambda\": %.17g,\n", params.lambda);
    std::fprintf(fp, "    \"T\": %.17g,\n", params.T);
    std::fprintf(fp, "    \"eta\": %.17g,\n", params.eta);
    std::fprintf(fp, "    \"Vg\": %.17g,\n", params.Vg);
    if (std::isinf(params.tau)) {
        std::fprintf(fp, "    \"tau\": \"Inf\"\n");
    } else {
        std::fprintf(fp, "    \"tau\": %.17g\n", params.tau);
    }
    std::fprintf(fp, "  },\n");
    std::fprintf(fp, "  \"bias_sweep\": {\n");
    std::fprintf(fp, "    \"Vsd_start\": %.17g,\n", params.Vsd_start);
    std::fprintf(fp, "    \"Vsd_end\": %.17g,\n", params.Vsd_end);
    std::fprintf(fp, "    \"Vsd_step\": %.17g\n", params.Vsd_step);
    std::fprintf(fp, "  },\n");

    auto write_array = [&](const char* name, const std::vector<double>& arr, bool comma) {
        std::fprintf(fp, "  \"%s\": [\n", name);
        for (int i = 0; i < n; ++i) {
            std::fprintf(fp, "    %.17g%s\n", arr[i], (i + 1 < n) ? "," : "");
        }
        std::fprintf(fp, "  ]%s\n", comma ? "," : "");
    };

    write_array("Vsd", res.Vsd, true);
    write_array("I_tol", res.I_tol, true);
    write_array("I_seq", res.I_seq, true);
    write_array("I_cot", res.I_cot, true);
    std::fprintf(fp, "  \"timestamp\": \"%s\"\n", tbuf);
    std::fprintf(fp, "}\n");

    std::fclose(fp);
}

void write_results_csv(const std::string& filepath, const SimulationResult& res) {
    FILE* fp = std::fopen(filepath.c_str(), "w");
    if (fp == nullptr) {
        return;
    }

    int n = static_cast<int>(res.Vsd.size());
    std::fprintf(fp, "Vsd_V,I_tol_A,I_seq_A,I_cot_A\n");
    for (int i = 0; i < n; ++i) {
        std::fprintf(fp, "%.6e,%.6e,%.6e,%.6e\n",
                     res.Vsd[i], res.I_tol[i], res.I_seq[i], res.I_cot[i]);
    }

    std::fclose(fp);
}

}
