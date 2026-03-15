use num_complex::Complex64;
use std::f64::consts::PI;

/// Even Bernoulli numbers B_2, B_4, ..., B_40.
static BERNOULLI_EVEN: [f64; 20] = [
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
    -261082718496449122051.0 / 13530.0,
];

/// Complex digamma function psi(z) via asymptotic series.
/// Reflection formula for Re(z) <= 0, recurrence shift until |z| >= 20,
/// then 20-term Bernoulli asymptotic expansion.
pub fn digamma_c(z: Complex64) -> Complex64 {
    let one = Complex64::new(1.0, 0.0);
    let half = Complex64::new(0.5, 0.0);

    let mut z_work = z;
    let mut result = Complex64::new(0.0, 0.0);
    let mut reflection = false;

    if z_work.re <= 0.0 {
        reflection = true;
        z_work = one - z_work;
    }

    while z_work.norm() < 20.0 {
        result -= one / z_work;
        z_work += one;
    }

    result += z_work.ln() - half / z_work;

    let z_sq = z_work * z_work;
    let mut power = z_sq;
    for k in 0..20 {
        result -= Complex64::from(BERNOULLI_EVEN[k] / (2.0 * (k + 1) as f64)) / power;
        power *= z_sq;
    }

    if reflection {
        let pi_z = Complex64::from(PI) * z;
        result -= Complex64::from(PI) * pi_z.cos() / pi_z.sin();
    }

    result
}

/// Complex trigamma function psi'(z) via asymptotic series.
/// Reflection formula for Re(z) <= 0, recurrence shift until |z| >= 10,
/// then 10-term Bernoulli asymptotic expansion.
pub fn trigamma_c(z: Complex64) -> Complex64 {
    let one = Complex64::new(1.0, 0.0);
    let half = Complex64::new(0.5, 0.0);

    let mut z_work = z;
    let mut result = Complex64::new(0.0, 0.0);
    let mut reflection = false;

    if z_work.re <= 0.0 {
        reflection = true;
        z_work = one - z_work;
    }

    while z_work.norm() < 10.0 {
        result += one / (z_work * z_work);
        z_work += one;
    }

    let iz = one / z_work;
    let iz2 = iz * iz;
    result += iz + half * iz2;

    let mut power = iz2 * iz;
    for k in 0..10 {
        result += Complex64::from(BERNOULLI_EVEN[k]) * power;
        power *= iz2;
    }

    if reflection {
        let sinval = (Complex64::from(PI) * z).sin();
        let ratio = Complex64::from(PI) / sinval;
        return ratio * ratio - result;
    }

    result
}
