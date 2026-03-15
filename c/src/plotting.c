#include "plotting.h"

#include <stdio.h>

void plot_iv(const char *csv_path, const char *pdf_path, const char *png_path,
             int N, double lambda, double T, double vmode,
             double alphaL, double alphaR, double eta, double Vg,
             double wall_time, const char *spec) {
    (void)spec;
    if (csv_path == NULL || pdf_path == NULL || png_path == NULL) {
        return;
    }

    FILE *gp = popen("gnuplot", "w");
    if (gp == NULL) {
        fprintf(stderr, "Warning: gnuplot not available, skipping plot generation.\n");
        return;
    }

    fprintf(gp, "set datafile separator ','\n");
    fprintf(gp, "set key left top box opaque\n");
    fprintf(gp, "set grid lw 0.5 lc rgb '#cccccc'\n");
    fprintf(gp, "set xlabel 'V_{sd} (V)'\n");
    fprintf(gp, "set ylabel '{/Italic I} (A)'\n");
    fprintf(gp,
            "set title 'N=%d, lambda=%.1f, T=%.1f K, hbar*omega=%.0f meV'\n",
            N, lambda, T, vmode * 1e3);
    fprintf(gp,
            "set label 1 'C: alpha_L=%.2f, alpha_R=%.2f, eta=%.1f, Vg=%.1f V\\n"
            "t=%.2f s' at graph 0.98,0.03 right front\n",
            alphaL, alphaR, eta, Vg, wall_time);
    fprintf(gp, "set style line 1 lc rgb '#111111' pt 7 ps 0.35 lw 1.0\n");
    fprintf(gp, "set style line 2 lc rgb '#0059b3' lw 1.5\n");
    fprintf(gp, "set style line 3 lc rgb '#cc1f1f' lw 1.5 dt 2\n");

    fprintf(gp, "set terminal pdfcairo enhanced font 'Times,10' size 12cm,9cm\n");
    fprintf(gp, "set output '%s'\n", pdf_path);
    fprintf(gp,
            "plot '%s' every ::1 using 1:2 with points ls 1 title 'I_{total}',"
            " '' every ::1 using 1:3 with lines ls 2 title 'I_{seq}',"
            " '' every ::1 using 1:4 with lines ls 3 title 'I_{cot}'\n",
            csv_path);

    fprintf(gp, "set terminal pngcairo enhanced font 'Times,10' size 12cm,9cm\n");
    fprintf(gp, "set output '%s'\n", png_path);
    fprintf(gp, "replot\n");
    fprintf(gp, "unset output\n");

    int status = pclose(gp);
    if (status != 0) {
        fprintf(stderr, "Warning: gnuplot failed to generate plots.\n");
    }
}
