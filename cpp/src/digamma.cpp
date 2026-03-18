#include "digamma.hpp"

#include <array>
#include <cmath>

namespace fc {

// B_2, B_4, ..., B_40 (even Bernoulli numbers)
static constexpr std::array<double, 20> BERNOULLI_EVEN = {
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

static constexpr double PI = 3.14159265358979323846;

Complex digamma(Complex z) {
    Complex z_work = z;
    Complex result(0.0, 0.0);
    bool reflection = false;

    if (z_work.real() <= 0.0) {
        reflection = true;
        z_work = 1.0 - z_work;
    }

    while (std::abs(z_work) < 20.0) {
        result -= 1.0 / z_work;
        z_work += 1.0;
    }

    result += std::log(z_work) - 1.0 / (2.0 * z_work);
    Complex z_sq = z_work * z_work;
    Complex power = z_sq;
    for (int k = 0; k < 20; ++k) {
        result -= BERNOULLI_EVEN[k] / (2.0 * static_cast<double>(k + 1) * power);
        power *= z_sq;
    }

    if (reflection) {
        return result - PI * std::cos(PI * z) / std::sin(PI * z);
    }
    return result;
}

Complex trigamma(Complex z) {
    Complex z_work = z;
    Complex result(0.0, 0.0);
    bool reflection = false;

    if (z_work.real() <= 0.0) {
        reflection = true;
        z_work = 1.0 - z_work;
    }

    while (std::abs(z_work) < 10.0) {
        result += 1.0 / (z_work * z_work);
        z_work += 1.0;
    }

    Complex iz = 1.0 / z_work;
    Complex iz2 = iz * iz;
    result += iz + 0.5 * iz2;
    Complex power = iz2 * iz;
    for (int k = 0; k < 10; ++k) {
        result += BERNOULLI_EVEN[k] * power;
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
