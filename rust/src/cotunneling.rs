use num_complex::Complex64;

use crate::fc_matrix::FCCache;
use crate::regularized::{regularized_i, regularized_j, regularized_j_matrix};

#[inline(always)]
fn sanitize(x: f64) -> f64 {
    if x.is_finite() { x } else { 0.0 }
}

fn rel_diff_log10(a: f64, b: f64) -> f64 {
    let z = Complex64::new(b, 0.0) / Complex64::new(a, 0.0);
    (z.ln() / 10.0_f64.ln()).norm()
}

fn rel_diff_log(a: f64, b: f64) -> f64 {
    let z = Complex64::new(b, 0.0) / Complex64::new(a, 0.0);
    z.ln().norm()
}

fn any_greater_than(arr: &[f64], threshold: f64) -> bool {
    arr.iter().any(|&x| x > threshold)
}

// ============================================================================
// Inner computation functions (fixed N)
// ============================================================================

/// Inner single sum for n=0→0 cotunneling.
pub fn m_sum_mmr(
    n: usize,
    q1: i32,
    q2_vec: &[i32],
    mu_l: f64,
    mu_r: f64,
    vmode: f64,
    epsilond: f64,
    t: f64,
    fc: &mut FCCache,
    out: &mut [f64],
) {
    let nq2 = q2_vec.len();
    if n == 0 || nq2 == 0 {
        return;
    }

    let mut mm_sq = vec![0.0; n * nq2];
    let e2_vec: Vec<f64> = q2_vec.iter().map(|&q2| mu_r - (q1 - q2) as f64 * vmode).collect();
    let eps_vec: Vec<f64> = (0..n as i32).map(|r| epsilond - (q1 - r) as f64 * vmode).collect();

    for r in 0..n {
        let fc_q1r = fc.get(q1, r as i32);
        for j in 0..nq2 {
            let fc_q2r = fc.get(q2_vec[j], r as i32);
            let prod = fc_q2r * fc_q1r;
            mm_sq[r * nq2 + j] = prod * prod;
        }
    }

    let mut jr = vec![0.0; n * nq2];
    regularized_j(mu_l, &e2_vec, &eps_vec, t, &mut jr);

    for j in 0..nq2 {
        let mut sum = 0.0;
        for r in 0..n {
            let term = mm_sq[r * nq2 + j] * jr[r * nq2 + j];
            sum += sanitize(term);
        }
        out[j] = sanitize(sum);
    }
}

/// Inner single sum for n=1→1 cotunneling.
pub fn m_sum_mmr11(
    n: usize,
    q1: i32,
    q2_vec: &[i32],
    mu_l: f64,
    mu_r: f64,
    vmode: f64,
    epsilond: f64,
    t: f64,
    fc: &mut FCCache,
    out: &mut [f64],
) {
    let nq2 = q2_vec.len();
    if n == 0 || nq2 == 0 {
        return;
    }

    let mut mm_sq = vec![0.0; n * nq2];
    let e2_vec: Vec<f64> = q2_vec.iter().map(|&q2| mu_r - (q1 - q2) as f64 * vmode).collect();
    let mut eps_matrix = vec![0.0; nq2 * n];

    for r in 0..n {
        let fc_q1r = fc.get(q1, r as i32);
        for j in 0..nq2 {
            let fc_q2r = fc.get(q2_vec[j], r as i32);
            let prod = fc_q2r * fc_q1r;
            mm_sq[r * nq2 + j] = prod * prod;
            eps_matrix[j * n + r] = epsilond + (q2_vec[j] - r as i32) as f64 * vmode;
        }
    }

    let mut jr = vec![0.0; n * nq2];
    regularized_j_matrix(mu_l, &e2_vec, &eps_matrix, nq2, n, t, &mut jr);

    for j in 0..nq2 {
        let mut sum = 0.0;
        for r in 0..n {
            let term = mm_sq[r * nq2 + j] * jr[r * nq2 + j];
            sum += sanitize(term);
        }
        out[j] = sanitize(sum);
    }
}

/// Inner double sum for n=0→0 cotunneling.
pub fn m_sum_mmmmrs(
    n: usize,
    q1: i32,
    q2: i32,
    mu_l: f64,
    mu_r: f64,
    vmode: f64,
    epsilond: f64,
    t: f64,
    fc: &mut FCCache,
) -> f64 {
    if n == 0 {
        return 0.0;
    }

    let eps1: Vec<f64> = (0..n as i32).map(|r| epsilond - (q1 - r) as f64 * vmode).collect();
    let eps2 = eps1.clone();

    let e2 = mu_r - (q1 - q2) as f64 * vmode;
    let mut irs = vec![0.0; n * n];
    regularized_i(mu_l, e2, &eps1, &eps2, t, &mut irs);

    let mut sum = 0.0;
    for r in 0..n {
        let fc_q2r = fc.get(q2, r as i32);
        let fc_q1r = fc.get(q1, r as i32);
        for s in 0..n {
            if r == s {
                continue;
            }
            let fc_q2s = fc.get(q2, s as i32);
            let fc_q1s = fc.get(q1, s as i32);
            let mmmm = fc_q2r * fc_q1r * fc_q2s * fc_q1s;
            let term = mmmm * irs[r * n + s];
            sum += sanitize(term);
        }
    }

    sanitize(sum)
}

/// Inner double sum for n=1→1 cotunneling.
pub fn m_sum_mmmmrs11(
    n: usize,
    q1: i32,
    q2: i32,
    mu_l: f64,
    mu_r: f64,
    vmode: f64,
    epsilond: f64,
    t: f64,
    fc: &mut FCCache,
) -> f64 {
    if n == 0 {
        return 0.0;
    }

    let eps1: Vec<f64> = (0..n as i32).map(|r| epsilond + (q2 - r) as f64 * vmode).collect();
    let eps2 = eps1.clone();

    let e2 = mu_r - (q1 - q2) as f64 * vmode;
    let mut irs = vec![0.0; n * n];
    regularized_i(mu_l, e2, &eps1, &eps2, t, &mut irs);

    let mut sum = 0.0;
    for r in 0..n {
        let fc_q2r = fc.get(q2, r as i32);
        let fc_q1r = fc.get(q1, r as i32);
        for s in 0..n {
            if r == s {
                continue;
            }
            let fc_q2s = fc.get(q2, s as i32);
            let fc_q1s = fc.get(q1, s as i32);
            let mmmm = fc_q2r * fc_q1r * fc_q2s * fc_q1s;
            let term = mmmm * irs[r * n + s];
            sum += sanitize(term);
        }
    }

    sanitize(sum)
}

// ============================================================================
// Convergence wrappers (adaptively increase N until converged)
// ============================================================================

/// Convergence wrapper for m_sum_mmr (n=0→0 single sum).
pub fn sum_mmr(
    q1: i32,
    q2_vec: &[i32],
    lambda: f64,
    mu_l: f64,
    mu_r: f64,
    vmode: f64,
    epsilond: f64,
    t: f64,
    fc: &mut FCCache,
    out: &mut [f64],
) {
    let nq2 = q2_vec.len();
    if nq2 == 0 {
        return;
    }

    let mut temp_n = (lambda.powf(2.2) * 3.0).round() as usize;
    let epsilon_conv = 1e-14;

    let mut temp_tol = vec![0.0; nq2];
    m_sum_mmr(temp_n, q1, q2_vec, mu_l, mu_r, vmode, epsilond, t, fc, &mut temp_tol);
    temp_n += 5;
    let mut temp_tol2 = vec![0.0; nq2];
    m_sum_mmr(temp_n, q1, q2_vec, mu_l, mu_r, vmode, epsilond, t, fc, &mut temp_tol2);

    let mut relative_diff: Vec<f64> = (0..nq2)
        .map(|j| rel_diff_log10(temp_tol[j], temp_tol2[j]))
        .collect();

    while any_greater_than(&relative_diff, epsilon_conv) {
        let loc_indices: Vec<usize> = (0..nq2)
            .filter(|&j| relative_diff[j] > epsilon_conv)
            .collect();
        if loc_indices.is_empty() {
            break;
        }

        let step = ((temp_n as f64 * 0.5).round() as usize).clamp(10, 20);
        temp_n += step;

        temp_tol.copy_from_slice(&temp_tol2);

        let q2_subset: Vec<i32> = loc_indices.iter().map(|&k| q2_vec[k]).collect();
        let mut partial = vec![0.0; loc_indices.len()];
        m_sum_mmr(temp_n, q1, &q2_subset, mu_l, mu_r, vmode, epsilond, t, fc, &mut partial);

        for (k, &idx) in loc_indices.iter().enumerate() {
            temp_tol2[idx] = partial[k];
        }

        relative_diff = (0..nq2)
            .map(|j| rel_diff_log10(temp_tol[j], temp_tol2[j]))
            .collect();
    }

    for j in 0..nq2 {
        out[j] = sanitize(temp_tol2[j]);
    }
}

/// Convergence wrapper for m_sum_mmr11 (n=1→1 single sum).
pub fn sum_mmr11(
    q1: i32,
    q2_vec: &[i32],
    lambda: f64,
    mu_l: f64,
    mu_r: f64,
    vmode: f64,
    epsilond: f64,
    t: f64,
    fc: &mut FCCache,
    out: &mut [f64],
) {
    let nq2 = q2_vec.len();
    if nq2 == 0 {
        return;
    }

    let mut temp_n = (lambda.powf(2.2) * 3.0).round() as usize;
    let epsilon_conv = 1e-14;

    let mut temp_tol = vec![0.0; nq2];
    m_sum_mmr11(temp_n, q1, q2_vec, mu_l, mu_r, vmode, epsilond, t, fc, &mut temp_tol);
    temp_n += 5;
    let mut temp_tol2 = vec![0.0; nq2];
    m_sum_mmr11(temp_n, q1, q2_vec, mu_l, mu_r, vmode, epsilond, t, fc, &mut temp_tol2);

    let mut relative_diff: Vec<f64> = (0..nq2)
        .map(|j| rel_diff_log10(temp_tol[j], temp_tol2[j]))
        .collect();

    while any_greater_than(&relative_diff, epsilon_conv) {
        let loc_indices: Vec<usize> = (0..nq2)
            .filter(|&j| relative_diff[j] > epsilon_conv)
            .collect();
        if loc_indices.is_empty() {
            break;
        }

        let step = ((temp_n as f64 * 0.5).round() as usize).clamp(10, 20);
        temp_n += step;

        temp_tol.copy_from_slice(&temp_tol2);

        let q2_subset: Vec<i32> = loc_indices.iter().map(|&k| q2_vec[k]).collect();
        let mut partial = vec![0.0; loc_indices.len()];
        m_sum_mmr11(temp_n, q1, &q2_subset, mu_l, mu_r, vmode, epsilond, t, fc, &mut partial);

        for (k, &idx) in loc_indices.iter().enumerate() {
            temp_tol2[idx] = partial[k];
        }

        relative_diff = (0..nq2)
            .map(|j| rel_diff_log10(temp_tol[j], temp_tol2[j]))
            .collect();
    }

    for j in 0..nq2 {
        out[j] = sanitize(temp_tol2[j]);
    }
}

/// Convergence wrapper for m_sum_mmmmrs (n=0→0 double sum).
pub fn sum_mmmmrs(
    q1: i32,
    q2_vec: &[i32],
    lambda: f64,
    mu_l: f64,
    mu_r: f64,
    vmode: f64,
    epsilond: f64,
    t: f64,
    fc: &mut FCCache,
    out: &mut [f64],
) {
    let nq2 = q2_vec.len();
    for idx in 0..nq2 {
        let q2 = q2_vec[idx];
        let mut temp_n = (lambda * lambda * 4.0).round() as usize;
        let epsilon_conv = 1e-14;

        let mut temp_tol = m_sum_mmmmrs(temp_n, q1, q2, mu_l, mu_r, vmode, epsilond, t, fc);
        temp_n += 5;
        let mut temp_tol2 = m_sum_mmmmrs(temp_n, q1, q2, mu_l, mu_r, vmode, epsilond, t, fc);

        let mut rd = rel_diff_log(temp_tol, temp_tol2);

        while rd.is_finite() && rd > epsilon_conv {
            temp_tol = temp_tol2;
            let step = ((temp_n as f64 * 0.5).round() as usize).clamp(20, 40);
            temp_n += step;
            temp_tol2 = m_sum_mmmmrs(temp_n, q1, q2, mu_l, mu_r, vmode, epsilond, t, fc);
            rd = rel_diff_log(temp_tol, temp_tol2);
        }

        out[idx] = sanitize(temp_tol2);
    }
}

/// Convergence wrapper for m_sum_mmmmrs11 (n=1→1 double sum).
pub fn sum_mmmmrs11(
    q1: i32,
    q2_vec: &[i32],
    lambda: f64,
    mu_l: f64,
    mu_r: f64,
    vmode: f64,
    epsilond: f64,
    t: f64,
    fc: &mut FCCache,
    out: &mut [f64],
) {
    let nq2 = q2_vec.len();
    for idx in 0..nq2 {
        let q2 = q2_vec[idx];
        let mut temp_n = (lambda * lambda * 4.0).round() as usize;
        let epsilon_conv = 1e-14;

        let mut temp_tol = m_sum_mmmmrs11(temp_n, q1, q2, mu_l, mu_r, vmode, epsilond, t, fc);
        temp_n += 5;
        let second_n = if temp_n > 1 { temp_n - 1 } else { 2 };
        let mut temp_tol2 =
            m_sum_mmmmrs11(second_n, q1, q2, mu_l, mu_r, vmode, epsilond, t, fc);

        let mut rd = rel_diff_log(temp_tol, temp_tol2);

        while rd.is_finite() && rd > epsilon_conv {
            let step = ((temp_n as f64 * 0.5).round() as usize).clamp(20, 40);
            temp_n += step;
            temp_tol = temp_tol2;
            temp_tol2 = m_sum_mmmmrs11(temp_n, q1, q2, mu_l, mu_r, vmode, epsilond, t, fc);
            rd = rel_diff_log(temp_tol, temp_tol2);
        }

        out[idx] = sanitize(temp_tol2);
    }
}
