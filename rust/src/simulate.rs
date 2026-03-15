use crate::constants::ELEMENTARY_CHARGE;
use crate::current::{current_from_rate_equations, CurrentResult};
use crate::fc_matrix::FCCache;
use crate::rate::{calculate_all_rate_w, RateStore};

pub struct SimulationResult {
    pub vsd: Vec<f64>,
    pub i_tol: Vec<f64>,
    pub i_seq: Vec<f64>,
    pub i_cot: Vec<f64>,
}

/// Run the full I-V simulation over all bias points.
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
    let mut result = SimulationResult {
        vsd: vsd_vec.to_vec(),
        i_tol: vec![0.0; n_vsd],
        i_seq: vec![0.0; n_vsd],
        i_cot: vec![0.0; n_vsd],
    };

    let mut fc = FCCache::new(lambda);

    for vv in 0..n_vsd {
        let v = vsd_vec[vv];
        if verbose {
            eprint!("\rBias point {}/{} (Vsd = {:.4} V)", vv + 1, n_vsd, v);
        }

        let mut store = RateStore::new(n);

        // Compute rates for both leads
        calculate_all_rate_w(
            &mut store, n, vmode, alpha_l, alpha_r, lambda, v, t, eta, 1, vg, &mut fc,
        );
        calculate_all_rate_w(
            &mut store, n, vmode, alpha_l, alpha_r, lambda, v, t, eta, -1, vg, &mut fc,
        );

        let mut cr: CurrentResult = current_from_rate_equations(n, vmode, t, tau, &store);

        // Sign correction: I *= -sign(Vsd)
        if v != 0.0 {
            let s = if v > 0.0 { -1.0 } else { 1.0 };
            cr.i_tol *= s;
            cr.i_seq *= s;
            cr.i_cot *= s;
        }

        // Convert to SI (Amperes)
        result.i_tol[vv] = cr.i_tol * ELEMENTARY_CHARGE;
        result.i_seq[vv] = cr.i_seq * ELEMENTARY_CHARGE;
        result.i_cot[vv] = cr.i_cot * ELEMENTARY_CHARGE;
    }

    if verbose {
        eprintln!();
    }

    result
}
