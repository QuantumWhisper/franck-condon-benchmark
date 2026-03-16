use num_complex::Complex64;

use crate::digamma_table::DigammaTable;
use crate::fc_matrix::FCCache;
use crate::regularized::{regularized_i, regularized_j, regularized_j_matrix};

const MAX_CONV_ITER: usize = 200;

#[inline(always)]
fn sanitize(x: f64) -> f64 {
    if x.is_finite() { x } else { 0.0 }
}

fn rel_diff_log10(a: f64, b: f64) -> f64 {
    if a.abs() < 1e-300 && b.abs() < 1e-300 {
        return 0.0;
    }
    if a.abs() < 1e-300 || b.abs() < 1e-300 {
        return f64::INFINITY;
    }
    let z = Complex64::new(b, 0.0) / Complex64::new(a, 0.0);
    (z.ln() / 10.0_f64.ln()).norm()
}

fn rel_diff_log(a: f64, b: f64) -> f64 {
    if a.abs() < 1e-300 && b.abs() < 1e-300 {
        return 0.0;
    }
    if a.abs() < 1e-300 || b.abs() < 1e-300 {
        return f64::INFINITY;
    }
    let z = Complex64::new(b, 0.0) / Complex64::new(a, 0.0);
    z.ln().norm()
}

fn any_greater_than(arr: &[f64], threshold: f64) -> bool {
    arr.iter().any(|&x| x > threshold)
}

pub fn m_sum_mmr(
    n: usize,
    q1: i32,
    q2_vec: &[i32],
    mu_l: f64,
    mu_r: f64,
    vmode: f64,
    epsilond: f64,
    t: f64,
    fc: &FCCache,
    out: &mut [f64],
    table: Option<&DigammaTable>,
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
    regularized_j(mu_l, &e2_vec, &eps_vec, t, &mut jr, table);

    for j in 0..nq2 {
        let mut sum = 0.0;
        for r in 0..n {
            let term = mm_sq[r * nq2 + j] * jr[r * nq2 + j];
            sum += sanitize(term);
        }
        out[j] = sanitize(sum);
    }
}

pub fn m_sum_mmr11(
    n: usize,
    q1: i32,
    q2_vec: &[i32],
    mu_l: f64,
    mu_r: f64,
    vmode: f64,
    epsilond: f64,
    t: f64,
    fc: &FCCache,
    out: &mut [f64],
    table: Option<&DigammaTable>,
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
    regularized_j_matrix(mu_l, &e2_vec, &eps_matrix, nq2, n, t, &mut jr, table);

    for j in 0..nq2 {
        let mut sum = 0.0;
        for r in 0..n {
            let term = mm_sq[r * nq2 + j] * jr[r * nq2 + j];
            sum += sanitize(term);
        }
        out[j] = sanitize(sum);
    }
}

pub fn m_sum_mmmmrs(
    n: usize,
    q1: i32,
    q2: i32,
    mu_l: f64,
    mu_r: f64,
    vmode: f64,
    epsilond: f64,
    t: f64,
    fc: &FCCache,
    table: Option<&DigammaTable>,
) -> f64 {
    if n == 0 {
        return 0.0;
    }

    let eps1: Vec<f64> = (0..n as i32).map(|r| epsilond - (q1 - r) as f64 * vmode).collect();
    let eps2 = eps1.clone();

    let e2 = mu_r - (q1 - q2) as f64 * vmode;
    let mut irs = vec![0.0; n * n];
    regularized_i(mu_l, e2, &eps1, &eps2, t, &mut irs, table);

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

pub fn m_sum_mmmmrs11(
    n: usize,
    q1: i32,
    q2: i32,
    mu_l: f64,
    mu_r: f64,
    vmode: f64,
    epsilond: f64,
    t: f64,
    fc: &FCCache,
    table: Option<&DigammaTable>,
) -> f64 {
    if n == 0 {
        return 0.0;
    }

    let eps1: Vec<f64> = (0..n as i32).map(|r| epsilond + (q2 - r) as f64 * vmode).collect();
    let eps2 = eps1.clone();

    let e2 = mu_r - (q1 - q2) as f64 * vmode;
    let mut irs = vec![0.0; n * n];
    regularized_i(mu_l, e2, &eps1, &eps2, t, &mut irs, table);

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

pub fn sum_mmr(
    q1: i32,
    q2_vec: &[i32],
    lambda: f64,
    mu_l: f64,
    mu_r: f64,
    vmode: f64,
    epsilond: f64,
    t: f64,
    fc: &FCCache,
    out: &mut [f64],
    table: Option<&DigammaTable>,
) {
    let nq2 = q2_vec.len();
    if nq2 == 0 {
        return;
    }

    let mut temp_n = (lambda.powf(2.2) * 3.0).round() as usize;
    let epsilon_conv = 1e-14;

    let mut temp_tol = vec![0.0; nq2];
    m_sum_mmr(temp_n, q1, q2_vec, mu_l, mu_r, vmode, epsilond, t, fc, &mut temp_tol, table);
    temp_n += 5;
    let mut temp_tol2 = vec![0.0; nq2];
    m_sum_mmr(temp_n, q1, q2_vec, mu_l, mu_r, vmode, epsilond, t, fc, &mut temp_tol2, table);

    let mut relative_diff: Vec<f64> = (0..nq2)
        .map(|j| rel_diff_log10(temp_tol[j], temp_tol2[j]))
        .collect();

    let mut conv_iter = 0usize;
    while any_greater_than(&relative_diff, epsilon_conv) && conv_iter < MAX_CONV_ITER {
        conv_iter += 1;
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
        m_sum_mmr(temp_n, q1, &q2_subset, mu_l, mu_r, vmode, epsilond, t, fc, &mut partial, table);

        for (k, &idx) in loc_indices.iter().enumerate() {
            temp_tol2[idx] = partial[k];
        }

        relative_diff = (0..nq2)
            .map(|j| rel_diff_log10(temp_tol[j], temp_tol2[j]))
            .collect();
    }
    if conv_iter >= MAX_CONV_ITER {
        eprintln!("WARNING: sum_mmr convergence not reached after {} iterations (q1={})", MAX_CONV_ITER, q1);
    }

    for j in 0..nq2 {
        out[j] = sanitize(temp_tol2[j]);
    }
}

pub fn sum_mmr11(
    q1: i32,
    q2_vec: &[i32],
    lambda: f64,
    mu_l: f64,
    mu_r: f64,
    vmode: f64,
    epsilond: f64,
    t: f64,
    fc: &FCCache,
    out: &mut [f64],
    table: Option<&DigammaTable>,
) {
    let nq2 = q2_vec.len();
    if nq2 == 0 {
        return;
    }

    let mut temp_n = (lambda.powf(2.2) * 3.0).round() as usize;
    let epsilon_conv = 1e-14;

    let mut temp_tol = vec![0.0; nq2];
    m_sum_mmr11(temp_n, q1, q2_vec, mu_l, mu_r, vmode, epsilond, t, fc, &mut temp_tol, table);
    temp_n += 5;
    let mut temp_tol2 = vec![0.0; nq2];
    m_sum_mmr11(temp_n, q1, q2_vec, mu_l, mu_r, vmode, epsilond, t, fc, &mut temp_tol2, table);

    let mut relative_diff: Vec<f64> = (0..nq2)
        .map(|j| rel_diff_log10(temp_tol[j], temp_tol2[j]))
        .collect();

    let mut conv_iter = 0usize;
    while any_greater_than(&relative_diff, epsilon_conv) && conv_iter < MAX_CONV_ITER {
        conv_iter += 1;
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
        m_sum_mmr11(temp_n, q1, &q2_subset, mu_l, mu_r, vmode, epsilond, t, fc, &mut partial, table);

        for (k, &idx) in loc_indices.iter().enumerate() {
            temp_tol2[idx] = partial[k];
        }

        relative_diff = (0..nq2)
            .map(|j| rel_diff_log10(temp_tol[j], temp_tol2[j]))
            .collect();
    }
    if conv_iter >= MAX_CONV_ITER {
        eprintln!("WARNING: sum_mmr11 convergence not reached after {} iterations (q1={})", MAX_CONV_ITER, q1);
    }

    for j in 0..nq2 {
        out[j] = sanitize(temp_tol2[j]);
    }
}

pub fn sum_mmmmrs(
    q1: i32,
    q2_vec: &[i32],
    lambda: f64,
    mu_l: f64,
    mu_r: f64,
    vmode: f64,
    epsilond: f64,
    t: f64,
    fc: &FCCache,
    out: &mut [f64],
    table: Option<&DigammaTable>,
) {
    let nq2 = q2_vec.len();
    for idx in 0..nq2 {
        let q2 = q2_vec[idx];
        let mut temp_n = (lambda * lambda * 4.0).round() as usize;
        let epsilon_conv = 1e-14;

        let mut temp_tol = m_sum_mmmmrs(temp_n, q1, q2, mu_l, mu_r, vmode, epsilond, t, fc, table);
        temp_n += 5;
        let mut temp_tol2 = m_sum_mmmmrs(temp_n, q1, q2, mu_l, mu_r, vmode, epsilond, t, fc, table);

        let mut rd = rel_diff_log(temp_tol, temp_tol2);
        let mut conv_iter = 0usize;

        while rd.is_finite() && rd > epsilon_conv && conv_iter < MAX_CONV_ITER {
            conv_iter += 1;
            temp_tol = temp_tol2;
            let step = ((temp_n as f64 * 0.5).round() as usize).clamp(20, 40);
            temp_n += step;
            temp_tol2 = m_sum_mmmmrs(temp_n, q1, q2, mu_l, mu_r, vmode, epsilond, t, fc, table);
            rd = rel_diff_log(temp_tol, temp_tol2);
        }
        if conv_iter >= MAX_CONV_ITER {
            eprintln!("WARNING: sum_mmmmrs convergence not reached after {} iterations (q1={}, q2={})", MAX_CONV_ITER, q1, q2);
        }

        out[idx] = sanitize(temp_tol2);
    }
}

pub fn sum_mmmmrs11(
    q1: i32,
    q2_vec: &[i32],
    lambda: f64,
    mu_l: f64,
    mu_r: f64,
    vmode: f64,
    epsilond: f64,
    t: f64,
    fc: &FCCache,
    out: &mut [f64],
    table: Option<&DigammaTable>,
) {
    let nq2 = q2_vec.len();
    for idx in 0..nq2 {
        let q2 = q2_vec[idx];
        let mut temp_n = (lambda * lambda * 4.0).round() as usize;
        let epsilon_conv = 1e-14;

        let mut temp_tol = m_sum_mmmmrs11(temp_n, q1, q2, mu_l, mu_r, vmode, epsilond, t, fc, table);
        temp_n += 5;
        let second_n = if temp_n > 1 { temp_n - 1 } else { 2 };
        let mut temp_tol2 =
            m_sum_mmmmrs11(second_n, q1, q2, mu_l, mu_r, vmode, epsilond, t, fc, table);

        let mut rd = rel_diff_log(temp_tol, temp_tol2);
        let mut conv_iter = 0usize;

        while rd.is_finite() && rd > epsilon_conv && conv_iter < MAX_CONV_ITER {
            conv_iter += 1;
            let step = ((temp_n as f64 * 0.5).round() as usize).clamp(20, 40);
            temp_n += step;
            temp_tol = temp_tol2;
            temp_tol2 = m_sum_mmmmrs11(temp_n, q1, q2, mu_l, mu_r, vmode, epsilond, t, fc, table);
            rd = rel_diff_log(temp_tol, temp_tol2);
        }
        if conv_iter >= MAX_CONV_ITER {
            eprintln!("WARNING: sum_mmmmrs11 convergence not reached after {} iterations (q1={}, q2={})", MAX_CONV_ITER, q1, q2);
        }

        out[idx] = sanitize(temp_tol2);
    }
}
