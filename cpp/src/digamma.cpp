#include "digamma.hpp"

#include <array>
#include <cmath>

namespace fc {

// B_2, B_4, ..., B_20 (10 terms). Used directly by trigamma asymptotic series.
static constexpr std::array<double, 10> BERNOULLI_EVEN = {
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

// DIGAMMA_COEFF[k] = B_{2(k+1)} / (2*(k+1)), pre-divided to eliminate runtime division
static constexpr std::array<double, 10> DIGAMMA_COEFF = {
     1.0 / 6.0 / 2.0,
    -1.0 / 30.0 / 4.0,
     1.0 / 42.0 / 6.0,
    -1.0 / 30.0 / 8.0,
     5.0 / 66.0 / 10.0,
    -691.0 / 2730.0 / 12.0,
     7.0 / 6.0 / 14.0,
    -3617.0 / 510.0 / 16.0,
     43867.0 / 798.0 / 18.0,
    -174611.0 / 330.0 / 20.0
};

static constexpr double PI = 3.14159265358979323846;

Complex digamma(Complex z) {
    // Fast path: 5-term series for |z|^2 > 900, Re(z) > 0. >99% hit rate at T=4.2K.
    if (z.real() > 0.0 && std::norm(z) > 900.0) {
        Complex inv_z = 1.0 / z;
        Complex inv_z_sq = inv_z * inv_z;
        Complex result = std::log(z) - inv_z * 0.5;
        Complex inv_power = inv_z_sq;
        for (int k = 0; k < 5; ++k) {
            result -= inv_power * DIGAMMA_COEFF[k];
            inv_power *= inv_z_sq;
        }
        return result;
    }

    Complex z_work = z;
    Complex result(0.0, 0.0);
    bool reflection = false;

    if (z_work.real() <= 0.0) {
        reflection = true;
        z_work = 1.0 - z_work;
    }

    // std::norm() returns |z|^2, avoiding sqrt per iteration
    while (std::norm(z_work) < 400.0) {
        result -= 1.0 / z_work;
        z_work += 1.0;
    }

    Complex inv_z = 1.0 / z_work;
    result += std::log(z_work) - inv_z * 0.5;
    Complex inv_z_sq = inv_z * inv_z;
    Complex inv_power = inv_z_sq;
    for (int k = 0; k < 10; ++k) {
        result -= inv_power * DIGAMMA_COEFF[k];
        inv_power *= inv_z_sq;
    }

    if (reflection) {
        return result - PI * std::cos(PI * z) / std::sin(PI * z);
    }
    return result;
}

Complex trigamma(Complex z) {
    // Fast path: 5-term series for |z|^2 > 900, Re(z) > 0
    if (z.real() > 0.0 && std::norm(z) > 900.0) {
        Complex iz = 1.0 / z;
        Complex iz2 = iz * iz;
        Complex result = iz + iz2 * 0.5;
        Complex power = iz2 * iz;
        for (int k = 0; k < 5; ++k) {
            result += power * BERNOULLI_EVEN[k];
            power *= iz2;
        }
        return result;
    }

    Complex z_work = z;
    Complex result(0.0, 0.0);
    bool reflection = false;

    if (z_work.real() <= 0.0) {
        reflection = true;
        z_work = 1.0 - z_work;
    }

    while (std::norm(z_work) < 100.0) {
        result += 1.0 / (z_work * z_work);
        z_work += 1.0;
    }

    Complex iz = 1.0 / z_work;
    Complex iz2 = iz * iz;
    result += iz + 0.5 * iz2;
    Complex power = iz2 * iz;
    for (int k = 0; k < 10; ++k) {
        result += power * BERNOULLI_EVEN[k];
        power *= iz2;
    }

    if (reflection) {
        Complex sinval = std::sin(PI * z);
        Complex ratio = PI / sinval;
        return ratio * ratio - result;
    }
    return result;
}

}
