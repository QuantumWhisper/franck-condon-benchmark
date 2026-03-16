use std::sync::atomic::{AtomicUsize, Ordering};

use rayon::prelude::*;

use crate::constants::ELEMENTARY_CHARGE;
use crate::current::current_from_rate_equations;
use crate::digamma_table::DigammaTable;
use crate::fc_matrix::{FCCache, FC_MAX_N};
use crate::rate::{calculate_all_rate_w, RateStore};

pub struct SimulationResult {
    pub vsd: Vec<f64>,
    pub i_tol: Vec<f64>,
    pub i_seq: Vec<f64>,
    pub i_cot: Vec<f64>,
}

pub fn simulate_iv(
    n: usize,
    vmode: f64,
    alpha_l: f64,
    alpha_r: f64,
    lambda: f64,
    vsd_vec: &[f64],
    t: f64,
    eta: f64,
    vg: f64,
    tau: f64,
    verbose: bool,
) -> SimulationResult {
    let n_vsd = vsd_vec.len();

    if n == 0 {
        return SimulationResult {
            vsd: vsd_vec.to_vec(),
            i_tol: vec![0.0; n_vsd],
            i_seq: vec![0.0; n_vsd],
            i_cot: vec![0.0; n_vsd],
        };
    }
    if t <= 0.0 {
        panic!("Temperature T must be positive (got T={:.6e}). T=0 is not supported.", t);
    }
    if vmode <= 0.0 {
        panic!("Phonon energy vmode must be positive (got vmode={:.6e}).", vmode);
    }
    if tau <= 0.0 && !tau.is_infinite() {
        panic!("Phonon relaxation time tau must be positive or Inf (got tau={:.6e}).", tau);
    }
    if alpha_l < 0.0 || alpha_r < 0.0 {
        panic!("Tunnel couplings alpha_L, alpha_R must be non-negative (got {:.6e}, {:.6e}).", alpha_l, alpha_r);
    }
    if lambda < 0.0 {
        panic!("Electron-phonon coupling lambda must be non-negative (got {:.6e}).", lambda);
    }
    let n_conv_est = (lambda.powf(2.2) * 3.0).round() as usize;
    if n_conv_est >= FC_MAX_N {
        eprintln!(
            "WARNING: lambda={:.1} requires convergence N~{} which exceeds FC cache size {}. \
             Results may be inaccurate. Consider reducing lambda or increasing FC_MAX_N.",
            lambda, n_conv_est, FC_MAX_N
        );
    }

    let mut fc_template = FCCache::new(lambda);
    let precompute_bound = ((lambda.powf(2.2) * 3.0).round() as usize + 50).max(n).min(FC_MAX_N);
    fc_template.precompute(precompute_bound);

    if verbose {
        eprint!("Building digamma lookup table...");
    }
    let dtable = DigammaTable::new();
    if verbose {
        eprintln!(" done.");
    }

    let progress = AtomicUsize::new(0);

    let results: Vec<(f64, f64, f64)> = vsd_vec
        .par_iter()
        .map(|&v| {
            let mut store = RateStore::new(n);

            calculate_all_rate_w(
                &mut store, n, vmode, alpha_l, alpha_r, lambda, v, t, eta, 1, vg, &fc_template,
                Some(&dtable),
            );
            calculate_all_rate_w(
                &mut store, n, vmode, alpha_l, alpha_r, lambda, v, t, eta, -1, vg, &fc_template,
                Some(&dtable),
            );

            let mut cr = current_from_rate_equations(n, vmode, t, tau, &store);

            if v != 0.0 {
                let s = if v > 0.0 { -1.0 } else { 1.0 };
                cr.i_tol *= s;
                cr.i_seq *= s;
                cr.i_cot *= s;
            }

            if verbose {
                let done = progress.fetch_add(1, Ordering::Relaxed) + 1;
                if done % 20 == 0 || done == n_vsd {
                    eprint!("\rBias point {}/{}", done, n_vsd);
                }
            }

            (
                cr.i_tol * ELEMENTARY_CHARGE,
                cr.i_seq * ELEMENTARY_CHARGE,
                cr.i_cot * ELEMENTARY_CHARGE,
            )
        })
        .collect();

    if verbose {
        eprintln!();
    }

    let mut result = SimulationResult {
        vsd: vsd_vec.to_vec(),
        i_tol: vec![0.0; n_vsd],
        i_seq: vec![0.0; n_vsd],
        i_cot: vec![0.0; n_vsd],
    };

    for (vv, &(i_tol, i_seq, i_cot)) in results.iter().enumerate() {
        result.i_tol[vv] = i_tol;
        result.i_seq[vv] = i_seq;
        result.i_cot[vv] = i_cot;
    }

    result
}
