#include "digamma.h"

#include <complex.h>
#include <math.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_result.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static const double BERNOULLI_EVEN[20] = {
    1.0 / 6.0,
    -1.0 / 30.0,
    1.0 / 42.0,
    -1.0 / 30.0,
    5.0 / 66.0,
    -691.0 / 2730.0,
    7.0 / 6.0,
    -3617.0 / 510.0,
    43867.0 / 798.0,
    -174611.0 / 330.0,
    854513.0 / 138.0,
    -236364091.0 / 2730.0,
    8553103.0 / 6.0,
    -23749461029.0 / 870.0,
    8615841276005.0 / 14322.0,
    -7709321041217.0 / 510.0,
    2577687858367.0 / 6.0,
    -26315271553053477373.0 / 1919190.0,
    2929993913841559.0 / 6.0,
    -261082718496449122051.0 / 13530.0
};

/* Use GSL's highly-optimized complex digamma implementation */
double _Complex digamma_c(double _Complex z) {
    gsl_sf_result re_result, im_result;
    int status = gsl_sf_complex_psi_e(creal(z), cimag(z), &re_result, &im_result);
    if (status != 0) {
        /* Fallback to our asymptotic series if GSL fails */
        double _Complex z_work = z;
        double _Complex result = 0.0 + 0.0 * I;
        int reflection = 0;
        if (creal(z_work) <= 0.0) {
            reflection = 1;
            z_work = 1.0 - z_work;
        }
        while (cabs(z_work) < 20.0) {
            result -= 1.0 / z_work;
            z_work += 1.0;
        }
        result += clog(z_work) - 1.0 / (2.0 * z_work);
        double _Complex z_sq = z_work * z_work;
        double _Complex power = z_sq;
        for (int k = 0; k < 20; ++k) {
            result -= BERNOULLI_EVEN[k] / (2.0 * (double)(k + 1) * power);
            power *= z_sq;
        }
        if (reflection) {
            return result - M_PI * ccos(M_PI * z) / csin(M_PI * z);
        }
        return result;
    }
    return re_result.val + I * im_result.val;
}

double _Complex trigamma_c(double _Complex z) {
    double _Complex z_work = z;
    double _Complex result = 0.0 + 0.0 * I;
    int reflection = 0;

    if (creal(z_work) <= 0.0) {
        reflection = 1;
        z_work = 1.0 - z_work;
    }

    while (cabs(z_work) < 10.0) {
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
