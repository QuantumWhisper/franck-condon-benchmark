use crate::constants::KB_EV;

#[inline]
pub fn fermi(x: f64, mu: f64, t: f64) -> f64 {
    1.0 / (((x - mu) / (KB_EV * t)).exp() + 1.0)
}

#[inline]
pub fn bose_fcn(x: f64, t: f64) -> f64 {
    let beta = 1.0 / (KB_EV * t);
    1.0 / ((x * beta).exp() - 1.0)
}
