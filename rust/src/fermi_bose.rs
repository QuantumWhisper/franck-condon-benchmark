use crate::constants::KB_EV;

#[inline]
pub fn fermi(x: f64, mu: f64, t: f64) -> f64 {
    if t <= 0.0 {
        return if x > mu {
            0.0
        } else if x < mu {
            1.0
        } else {
            0.5
        };
    }
    let arg = (x - mu) / (KB_EV * t);
    if arg > 500.0 {
        0.0
    } else if arg < -500.0 {
        1.0
    } else {
        1.0 / (arg.exp() + 1.0)
    }
}

#[inline]
pub fn bose_fcn(x: f64, t: f64) -> f64 {
    if t <= 0.0 {
        return if x > 0.0 { 0.0 } else { 0.0 };
    }
    let beta = 1.0 / (KB_EV * t);
    let arg = x * beta;
    if arg.abs() < 1e-15 {
        return beta / x.signum().max(1.0);
    }
    if arg > 500.0 {
        0.0
    } else {
        1.0 / (arg.exp() - 1.0)
    }
}
