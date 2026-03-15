use crate::matrix::generate_matrix_w;
use crate::rate::RateStore;
use crate::solver::solve_steady_state;

pub struct CurrentResult {
    pub i_tol: f64,
    pub i_seq: f64,
    pub i_cot: f64,
}

/// Compute sequential and cotunneling currents from steady-state probabilities.
pub fn current_from_rate_equations(
    n: usize,
    vmode: f64,
    t: f64,
    tau: f64,
    store: &RateStore,
) -> CurrentResult {
    let dim = 2 * n;

    let mut w = vec![0.0; dim * dim];
    generate_matrix_w(store, n, vmode, t, tau, &mut w);

    let p = solve_steady_state(&w, n);

    // Sequential current: n=0→1 transitions
    // Sign convention: diffW = w_R - w_L
    let mut i_seq0 = 0.0;
    for q1 in 0..n {
        let mut sum_diff = 0.0;
        for q2 in 0..n {
            let w_r = store.rate_w(0, 1, q1, q2, -1); // right lead
            let w_l = store.rate_w(0, 1, q1, q2, 1); // left lead
            sum_diff += w_r - w_l;
        }
        i_seq0 += p[q1] * sum_diff;
    }

    // Sequential current: n=1→0 transitions
    // Sign convention: diffW = w_L - w_R (OPPOSITE!)
    let mut i_seq1 = 0.0;
    for q1 in 0..n {
        let mut sum_diff = 0.0;
        for q2 in 0..n {
            let w_r = store.rate_w(1, 0, q1, q2, -1);
            let w_l = store.rate_w(1, 0, q1, q2, 1);
            sum_diff += w_l - w_r;
        }
        i_seq1 += p[n + q1] * sum_diff;
    }

    let i_seq = i_seq0 + i_seq1;

    // Cotunneling current for both n=0→0 and n=1→1
    // Sign convention: diffW = w_RL - w_LR
    let mut i_cot = 0.0;
    for nn in 0..=1_usize {
        for q1 in 0..n {
            let mut sum_diff = 0.0;
            for q2 in 0..n {
                let w_rl = store.rate_w(nn, nn, q1, q2, -1); // right lead
                let w_lr = store.rate_w(nn, nn, q1, q2, 1); // left lead
                sum_diff += w_rl - w_lr;
            }
            i_cot += p[nn * n + q1] * sum_diff;
        }
    }

    CurrentResult {
        i_tol: i_seq + i_cot,
        i_seq,
        i_cot,
    }
}
