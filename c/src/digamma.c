#include "digamma.h"

#include <complex.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Even Bernoulli numbers B_2, B_4, ..., B_20 (10 terms).
 * Used directly by trigamma asymptotic series. */
static const double BERNOULLI_EVEN[10] = {
    1.0 / 6.0,
    -1.0 / 30.0,
    1.0 / 42.0,
    -1.0 / 30.0,
    5.0 / 66.0,
    -691.0 / 2730.0,
    7.0 / 6.0,
    -3617.0 / 510.0,
    43867.0 / 798.0,
    -174611.0 / 330.0
};

/* Pre-computed digamma coefficients: B_{2(k+1)} / (2*(k+1)) for k=0..9.
 * Eliminates runtime division in the asymptotic series:
 *   psi(z) ~ ln(z) - 1/(2z) - sum_k coeff[k] * z^{-2(k+1)}
 * Matches rust/src/digamma.rs DIGAMMA_COEFF exactly. */
static const double DIGAMMA_COEFF[10] = {
    1.0 / 6.0 / 2.0,             /* B_2  / 2  */
    -1.0 / 30.0 / 4.0,           /* B_4  / 4  */
    1.0 / 42.0 / 6.0,            /* B_6  / 6  */
    -1.0 / 30.0 / 8.0,           /* B_8  / 8  */
    5.0 / 66.0 / 10.0,           /* B_10 / 10 */
    -691.0 / 2730.0 / 12.0,      /* B_12 / 12 */
    7.0 / 6.0 / 14.0,            /* B_14 / 14 */
    -3617.0 / 510.0 / 16.0,      /* B_16 / 16 */
    43867.0 / 798.0 / 18.0,      /* B_18 / 18 */
    -174611.0 / 330.0 / 20.0     /* B_20 / 20 */
};

/* Fast-path digamma for |z|^2 > 900 and Re(z) > 0.
 * Uses only 5 Bernoulli terms — no reflection or recurrence needed.
 * At T=4.2K, >99% of digamma calls hit this path.
 * See rust/src/digamma.rs digamma_asymptotic5(). */
static inline double _Complex digamma_asymptotic5(double _Complex z) {
    double _Complex inv_z = 1.0 / z;
    double _Complex inv_z_sq = inv_z * inv_z;
    double _Complex result = clog(z) - inv_z * 0.5;
    double _Complex inv_power = inv_z_sq;
    for (int k = 0; k < 5; ++k) {
        result -= inv_power * DIGAMMA_COEFF[k];
        inv_power *= inv_z_sq;
    }
    return result;
}

/* Fast-path trigamma for |z|^2 > 900 and Re(z) > 0.
 * Uses only 5 Bernoulli terms — no reflection or recurrence needed.
 * See rust/src/digamma.rs trigamma_asymptotic5(). */
static inline double _Complex trigamma_asymptotic5(double _Complex z) {
    double _Complex iz = 1.0 / z;
    double _Complex iz2 = iz * iz;
    double _Complex result = iz + iz2 * 0.5;
    double _Complex power = iz2 * iz;
    for (int k = 0; k < 5; ++k) {
        result += power * BERNOULLI_EVEN[k];
        power *= iz2;
    }
    return result;
}

/* Complex digamma function psi(z) via pure asymptotic series.
 * 10-term Bernoulli expansion with pre-computed coefficients.
 * Uses norm_sqr thresholds to avoid sqrt (cabs) per iteration.
 * Matches Rust/Fortran/C++ implementations exactly. */
double _Complex digamma_c(double _Complex z) {
    /* Fast path: Re(z) > 0 and |z|^2 > 900 — skip reflection and recurrence */
    double re = creal(z), im = cimag(z);
    double norm_sq = re * re + im * im;
    if (re > 0.0 && norm_sq > 900.0) {
        return digamma_asymptotic5(z);
    }

    double _Complex z_work = z;
    double _Complex result = 0.0 + 0.0 * I;
    int reflection = 0;

    if (creal(z_work) <= 0.0) {
        reflection = 1;
        z_work = 1.0 - z_work;
    }

    /* Recurrence shift until |z|^2 >= 100 (i.e., |z| >= 10) */
    while (creal(z_work) * creal(z_work) + cimag(z_work) * cimag(z_work) < 100.0) {
        result -= 1.0 / z_work;
        z_work += 1.0;
    }

    /* Asymptotic expansion: multiply-instead-of-divide pattern */
    double _Complex inv_z = 1.0 / z_work;
    double _Complex inv_z_sq = inv_z * inv_z;
    result += clog(z_work) - inv_z * 0.5;
    double _Complex inv_power = inv_z_sq;
    for (int k = 0; k < 10; ++k) {
        result -= inv_power * DIGAMMA_COEFF[k];
        inv_power *= inv_z_sq;
    }

    if (reflection) {
        return result - M_PI * ccos(M_PI * z) / csin(M_PI * z);
    }
    return result;
}

/* Complex trigamma function psi'(z) via asymptotic series.
 * 10-term Bernoulli expansion, threshold |z| >= 10.
 * Uses norm_sqr to avoid sqrt per recurrence iteration. */
double _Complex trigamma_c(double _Complex z) {
    /* Fast path: Re(z) > 0 and |z|^2 > 900 — skip reflection and recurrence */
    double re = creal(z), im = cimag(z);
    double norm_sq = re * re + im * im;
    if (re > 0.0 && norm_sq > 900.0) {
        return trigamma_asymptotic5(z);
    }

    double _Complex z_work = z;
    double _Complex result = 0.0 + 0.0 * I;
    int reflection = 0;

    if (creal(z_work) <= 0.0) {
        reflection = 1;
        z_work = 1.0 - z_work;
    }

    /* Recurrence shift until |z|^2 >= 100 (i.e., |z| >= 10) */
    while (creal(z_work) * creal(z_work) + cimag(z_work) * cimag(z_work) < 100.0) {
        result += 1.0 / (z_work * z_work);
        z_work += 1.0;
    }

    double _Complex iz = 1.0 / z_work;
    double _Complex iz2 = iz * iz;
    result += iz + 0.5 * iz2;
    double _Complex power = iz2 * iz;
    for (int k = 0; k < 10; ++k) {
        result += BERNOULLI_EVEN[k] * power;
        power *= iz2;
    }

    if (reflection) {
        double _Complex sinval = csin(M_PI * z);
        double _Complex ratio = M_PI / sinval;
        return ratio * ratio - result;
    }
    return result;
}
