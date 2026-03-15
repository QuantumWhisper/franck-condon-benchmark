#include "laguerre.h"

double laguerre_L(int n, int alpha, double x) {
    if (n <= 0) {
        return 1.0;
    }

    double L_prev = 1.0;
    double L_curr = 1.0 + (double)alpha - x;
    if (n == 1) {
        return L_curr;
    }

    for (int k = 1; k < n; ++k) {
        double L_next =
            ((2.0 * (double)k + 1.0 + (double)alpha - x) * L_curr -
             ((double)k + (double)alpha) * L_prev) /
            ((double)k + 1.0);
        L_prev = L_curr;
        L_curr = L_next;
    }

    return L_curr;
}
