#pragma once

#include <string>

namespace fc {

void plot_iv(const std::string& csv_path, const std::string& pdf_path,
             const std::string& png_path,
             int N, double lambda, double T, double vmode,
             double alphaL, double alphaR, double eta, double Vg,
             double wall_time, const std::string& spec);

}
