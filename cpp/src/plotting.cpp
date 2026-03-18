#include "plotting.hpp"

#include <cstdio>

namespace fc {

void plot_iv(const std::string& csv_path, const std::string& pdf_path,
             const std::string& png_path,
             int N, double lambda, double T, double vmode,
             double alphaL, double alphaR, double eta, double Vg,
             double wall_time, const std::string& spec) {
    (void)spec;

    FILE* gp = popen("gnuplot", "w");
    if (gp == nullptr) {
        std::fprintf(stderr, "Warning: gnuplot not available, skipping plot generation.\n");
        return;
    }

    std::fprintf(gp, "set datafile separator ','\n");
    std::fprintf(gp, "set key left top box opaque\n");
    std::fprintf(gp, "set grid lw 0.5 lc rgb '#cccccc'\n");
    std::fprintf(gp, "set xlabel 'V_{sd} (V)'\n");
    std::fprintf(gp, "set ylabel '{/Italic I} (A)'\n");
    std::fprintf(gp,
            "set title 'N=%d, lambda=%.1f, T=%.1f K, hbar*omega=%.0f meV'\n",
            N, lambda, T, vmode * 1e3);
    std::fprintf(gp,
            "set label 1 'C++: alpha_L=%.2f, alpha_R=%.2f, eta=%.1f, Vg=%.1f V\\n"
            "t=%.2f s' at graph 0.98,0.03 right front\n",
            alphaL, alphaR, eta, Vg, wall_time);
    std::fprintf(gp, "set style line 1 lc rgb '#111111' pt 7 ps 0.35 lw 1.0\n");
    std::fprintf(gp, "set style line 2 lc rgb '#0059b3' lw 1.5\n");
    std::fprintf(gp, "set style line 3 lc rgb '#cc1f1f' lw 1.5 dt 2\n");

    std::fprintf(gp, "set terminal pdfcairo enhanced font 'Times,10' size 12cm,9cm\n");
    std::fprintf(gp, "set output '%s'\n", pdf_path.c_str());
    std::fprintf(gp,
            "plot '%s' every ::1 using 1:2 with points ls 1 title 'I_{total}',"
            " '' every ::1 using 1:3 with lines ls 2 title 'I_{seq}',"
            " '' every ::1 using 1:4 with lines ls 3 title 'I_{cot}'\n",
            csv_path.c_str());

    std::fprintf(gp, "set terminal pngcairo enhanced font 'Times,10' size 12cm,9cm\n");
    std::fprintf(gp, "set output '%s'\n", png_path.c_str());
    std::fprintf(gp, "replot\n");
    std::fprintf(gp, "unset output\n");

    int status = pclose(gp);
    if (status != 0) {
        std::fprintf(stderr, "Warning: gnuplot failed to generate plots.\n");
    }
}

}
