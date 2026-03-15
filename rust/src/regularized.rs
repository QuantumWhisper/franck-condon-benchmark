use num_complex::Complex64;
use std::f64::consts::PI;

use crate::constants::KB_EV;
use crate::digamma::{digamma_c, trigamma_c};
use crate::fermi_bose::bose_fcn;

/// Regularized I integral with factored digamma precomputation.
///
/// Key optimization from C port: digamma arguments factor into row-only
/// and column-only terms. Precompute per row/col → 4*N calls instead of 4*N*N.
///
/// Output is stored in row-major order: out[i * n_eps2 + j].
pub fn regularized_i(
    e1: f64,
    e2: f64,
    epsilon1: &[f64],
    epsilon2: &[f64],
    t: f64,
    out: &mut [f64],
) {
    let n_eps1 = epsilon1.len();
    let n_eps2 = epsilon2.len();
    if n_eps1 == 0 || n_eps2 == 0 {
        return;
    }

    let beta = 1.0 / (KB_EV * t);
    let bose_val = bose_fcn(e2 - e1, t);
    let coeff = beta / (2.0 * PI);

    // Precompute row-only digamma differences: psi(a1) - psi(a3)
    let row_diff: Vec<Complex64> = epsilon1
        .iter()
        .map(|&eps| {
            let a1 = Complex64::new(0.5, coeff * (e2 - eps));
            let a3 = Complex64::new(0.5, coeff * (e1 - eps));
            digamma_c(a1) - digamma_c(a3)
        })
        .collect();

    // Precompute column-only digamma differences: psi(a2) - psi(a4)
    let col_diff: Vec<Complex64> = epsilon2
        .iter()
        .map(|&eps| {
            let a2 = Complex64::new(0.5, -coeff * (e2 - eps));
            let a4 = Complex64::new(0.5, -coeff * (e1 - eps));
            digamma_c(a2) - digamma_c(a4)
        })
        .collect();

    // Compute matrix entries
    for i in 0..n_eps1 {
        for j in 0..n_eps2 {
            let denom = epsilon1[i] - epsilon2[j];
            let val = bose_val / denom * (row_diff[i] - col_diff[j]).re;
            out[i * n_eps2 + j] = if val.is_finite() { val } else { 0.0 };
        }
    }
}

/// Regularized J integral (vector E2).
/// Output is stored in row-major order: out[i * n_e2 + j].
pub fn regularized_j(
    e1: f64,
    e2: &[f64],
    epsilon: &[f64],
    t: f64,
    out: &mut [f64],
) {
    let n_e2 = e2.len();
    let n_eps = epsilon.len();
    if n_e2 == 0 || n_eps == 0 {
        return;
    }

    let beta = 1.0 / (KB_EV * t);
    let coeff = beta / (2.0 * PI);

    let bose_vals: Vec<f64> = e2.iter().map(|&e| bose_fcn(e - e1, t)).collect();

    // Precompute trigamma(a2) per row
    let trig_a2: Vec<Complex64> = epsilon
        .iter()
        .map(|&eps| {
            let a2 = Complex64::new(0.5, coeff * (e1 - eps));
            trigamma_c(a2)
        })
        .collect();

    for i in 0..n_eps {
        for j in 0..n_e2 {
            let a1 = Complex64::new(0.5, coeff * (e2[j] - epsilon[i]));
            let val = coeff * bose_vals[j] * (trigamma_c(a1) - trig_a2[i]).im;
            out[i * n_e2 + j] = if val.is_finite() { val } else { 0.0 };
        }
    }
}

/// Regularized J integral (matrix epsilon) for n=1→1 cotunneling.
/// epsilon_matrix is stored in row-major: epsilon_matrix[j * eps_cols + i].
/// Output is stored: out[i * n_e2 + j].
pub fn regularized_j_matrix(
    e1: f64,
    e2: &[f64],
    epsilon_matrix: &[f64],
    eps_rows: usize,
    eps_cols: usize,
    t: f64,
    out: &mut [f64],
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
            let a1 = Complex64::new(0.5, coeff * (e2[j] - eps));
            let a2 = Complex64::new(0.5, coeff * (e1 - eps));

            let val = coeff * bose_vals[j] * (trigamma_c(a1) - trigamma_c(a2)).im;
            out[i * n_e2 + j] = if val.is_finite() { val } else { 0.0 };
        }
    }
}
