#ifndef FC_PLOTTING_H
#define FC_PLOTTING_H

void plot_iv(const char *csv_path, const char *pdf_path, const char *png_path,
             int N, double lambda, double T, double vmode,
             double alphaL, double alphaR, double eta, double Vg,
             double wall_time, const char *spec);

#endif
