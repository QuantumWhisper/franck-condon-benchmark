use num_complex::Complex64;
use std::f64::consts::PI;

/// Even Bernoulli numbers B_2, B_4, ..., B_20 (10 terms).
/// Used directly by trigamma asymptotic series.
static BERNOULLI_EVEN: [f64; 10] = [
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
];

/// Pre-computed digamma coefficients: B_{2(k+1)} / (2*(k+1)) for k=0..9.
/// The asymptotic series is: psi(z) ~ ln(z) - 1/(2z) - sum_{k} coeff[k] / z^{2(k+1)}
static DIGAMMA_COEFF: [f64; 10] = [
    1.0 / 6.0 / 2.0,             // B_2  / 2
    -1.0 / 30.0 / 4.0,           // B_4  / 4
    1.0 / 42.0 / 6.0,            // B_6  / 6
    -1.0 / 30.0 / 8.0,           // B_8  / 8
    5.0 / 66.0 / 10.0,           // B_10 / 10
    -691.0 / 2730.0 / 12.0,      // B_12 / 12
    7.0 / 6.0 / 14.0,            // B_14 / 14
    -3617.0 / 510.0 / 16.0,      // B_16 / 16
    43867.0 / 798.0 / 18.0,      // B_18 / 18
    -174611.0 / 330.0 / 20.0,    // B_20 / 20
];

#[inline(always)]
pub fn digamma_asymptotic5(z: Complex64) -> Complex64 {
    let inv_z = z.inv();
    let inv_z_sq = inv_z * inv_z;
    let mut result = z.ln() - inv_z * 0.5;
    let mut inv_power = inv_z_sq;
    for k in 0..5 {
        result -= inv_power * DIGAMMA_COEFF[k];
        inv_power *= inv_z_sq;
    }
    result
}

#[inline(always)]
pub fn trigamma_asymptotic5(z: Complex64) -> Complex64 {
    let iz = z.inv();
    let iz2 = iz * iz;
    let mut result = iz + iz2 * 0.5;
    let mut power = iz2 * iz;
    for k in 0..5 {
        result += power * BERNOULLI_EVEN[k];
        power *= iz2;
    }
    result
}

#[inline]
pub fn digamma_c(z: Complex64) -> Complex64 {
    let mut z_work = z;
    let mut result = Complex64::new(0.0, 0.0);
    let mut reflection = false;

    if z_work.re <= 0.0 {
        reflection = true;
        z_work = Complex64::new(1.0 - z_work.re, -z_work.im);
    }

    // Recurrence shift until |z|^2 >= 100 (i.e., |z| >= 10)
    while z_work.norm_sqr() < 100.0 {
        result -= z_work.inv();
        z_work += 1.0;
    }

    // Asymptotic expansion: psi(z) ~ ln(z) - 1/(2z) - sum coeff[k] * z^{-2(k+1)}
    let inv_z = z_work.inv();
    result += z_work.ln() - inv_z * 0.5;

    let inv_z_sq = inv_z * inv_z;
    let mut inv_power = inv_z_sq;
    for k in 0..10 {
        result -= inv_power * DIGAMMA_COEFF[k];
        inv_power *= inv_z_sq;
    }

    if reflection {
        let pi_z = z * PI;
        result -= pi_z.cos() / pi_z.sin() * PI;
    }

    result
}

/// Complex trigamma function psi'(z) via asymptotic series.
///
/// 10-term Bernoulli expansion, threshold |z| >= 10.
#[inline]
pub fn trigamma_c(z: Complex64) -> Complex64 {
    let mut z_work = z;
    let mut result = Complex64::new(0.0, 0.0);
    let mut reflection = false;

    if z_work.re <= 0.0 {
        reflection = true;
        z_work = Complex64::new(1.0 - z_work.re, -z_work.im);
    }

    // Recurrence shift until |z|^2 >= 100
    while z_work.norm_sqr() < 100.0 {
        result += (z_work * z_work).inv();
        z_work += 1.0;
    }

    let iz = z_work.inv();
    let iz2 = iz * iz;
    result += iz + iz2 * 0.5;

    let mut power = iz2 * iz;
    for k in 0..10 {
        result += power * BERNOULLI_EVEN[k];
        power *= iz2;
    }

    if reflection {
        let sinval = (z * PI).sin();
        let ratio = sinval.inv() * PI;
        return ratio * ratio - result;
    }

    result
}
