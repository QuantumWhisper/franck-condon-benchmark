use num_complex::Complex64;
use std::f64::consts::PI;

use crate::constants::KB_EV;
use crate::digamma::{digamma_asymptotic5, digamma_c, trigamma_asymptotic5, trigamma_c};
use crate::digamma_table::DigammaTable;
use crate::fermi_bose::bose_fcn;

#[inline(always)]
fn fast_digamma_half(y: f64) -> Complex64 {
    let z = Complex64::new(0.5, y);
    if y * y > 900.0 {
        digamma_asymptotic5(z)
    } else {
        digamma_c(z)
    }
}

#[inline(always)]
fn fast_trigamma_half(y: f64) -> Complex64 {
    let z = Complex64::new(0.5, y);
    if y * y > 900.0 {
        trigamma_asymptotic5(z)
    } else {
        trigamma_c(z)
    }
}

pub fn regularized_i(
    e1: f64,
    e2: f64,
    epsilon1: &[f64],
    epsilon2: &[f64],
    t: f64,
    out: &mut [f64],
    table: Option<&DigammaTable>,
) {
    let n_eps1 = epsilon1.len();
    let n_eps2 = epsilon2.len();
    if n_eps1 == 0 || n_eps2 == 0 {
        return;
    }

    let beta = 1.0 / (KB_EV * t);
    let bose_val = bose_fcn(e2 - e1, t);
    let coeff = beta / (2.0 * PI);

    let row_diff: Vec<Complex64> = epsilon1
        .iter()
        .map(|&eps| {
            let y1 = coeff * (e2 - eps);
            let y3 = coeff * (e1 - eps);
            if let Some(tbl) = table {
                tbl.digamma(y1) - tbl.digamma(y3)
            } else {
                fast_digamma_half(y1) - fast_digamma_half(y3)
            }
        })
        .collect();

    let col_diff: Vec<Complex64> = epsilon2
        .iter()
        .map(|&eps| {
            let y2 = -coeff * (e2 - eps);
            let y4 = -coeff * (e1 - eps);
            if let Some(tbl) = table {
                tbl.digamma(y2) - tbl.digamma(y4)
            } else {
                fast_digamma_half(y2) - fast_digamma_half(y4)
            }
        })
        .collect();

    for i in 0..n_eps1 {
        for j in 0..n_eps2 {
            let denom = epsilon1[i] - epsilon2[j];
            let val = bose_val / denom * (row_diff[i] - col_diff[j]).re;
            out[i * n_eps2 + j] = if val.is_finite() { val } else { 0.0 };
        }
    }
}

pub fn regularized_j(
    e1: f64,
    e2: &[f64],
    epsilon: &[f64],
    t: f64,
    out: &mut [f64],
    table: Option<&DigammaTable>,
) {
    let n_e2 = e2.len();
    let n_eps = epsilon.len();
    if n_e2 == 0 || n_eps == 0 {
        return;
    }

    let beta = 1.0 / (KB_EV * t);
    let coeff = beta / (2.0 * PI);

    let bose_vals: Vec<f64> = e2.iter().map(|&e| bose_fcn(e - e1, t)).collect();

    let trig_a2: Vec<Complex64> = epsilon
        .iter()
        .map(|&eps| {
            let y = coeff * (e1 - eps);
            if let Some(tbl) = table {
                tbl.trigamma(y)
            } else {
                fast_trigamma_half(y)
            }
        })
        .collect();

    for i in 0..n_eps {
        for j in 0..n_e2 {
            let y = coeff * (e2[j] - epsilon[i]);
            let trig_a1 = if let Some(tbl) = table {
                tbl.trigamma(y)
            } else {
                fast_trigamma_half(y)
            };
            let val = coeff * bose_vals[j] * (trig_a1 - trig_a2[i]).im;
            out[i * n_e2 + j] = if val.is_finite() { val } else { 0.0 };
        }
    }
}

pub fn regularized_j_matrix(
    e1: f64,
    e2: &[f64],
    epsilon_matrix: &[f64],
    eps_rows: usize,
    eps_cols: usize,
    t: f64,
    out: &mut [f64],
    table: Option<&DigammaTable>,
) {
    let n_e2 = e2.len();
    if n_e2 == 0 || eps_rows == 0 || eps_cols == 0 {
        return;
    }

    let beta = 1.0 / (KB_EV * t);
    let coeff = beta / (2.0 * PI);

    let bose_vals: Vec<f64> = e2.iter().map(|&e| bose_fcn(e - e1, t)).collect();

    for i in 0..eps_cols {
        for j in 0..n_e2 {
            if j >= eps_rows {
                out[i * n_e2 + j] = 0.0;
                continue;
            }

            let eps = epsilon_matrix[j * eps_cols + i];
            let y1 = coeff * (e2[j] - eps);
            let y2 = coeff * (e1 - eps);

            let (ta1, ta2) = if let Some(tbl) = table {
                (tbl.trigamma(y1), tbl.trigamma(y2))
            } else {
                (fast_trigamma_half(y1), fast_trigamma_half(y2))
            };

            let val = coeff * bose_vals[j] * (ta1 - ta2).im;
            out[i * n_e2 + j] = if val.is_finite() { val } else { 0.0 };
        }
    }
}
