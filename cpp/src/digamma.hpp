#pragma once

#include <complex>

namespace fc {

using Complex = std::complex<double>;

Complex digamma(Complex z);
Complex trigamma(Complex z);

}
