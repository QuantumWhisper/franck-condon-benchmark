use std::f64::consts::PI;

use crate::constants::HBAR_EV;
use crate::cotunneling::{sum_mmmmrs, sum_mmmmrs11, sum_mmr, sum_mmr11};
use crate::fc_matrix::FCCache;
use crate::fermi_bose::fermi;

/// Flat-array rate store indexed by [n1][n2][q1][lead_idx][q2].
pub struct RateStore {
    pub n: usize,
    pub data: Vec<f64>,
}

impl RateStore {
    pub fn new(n: usize) -> Self {
        RateStore {
            n,
            data: vec![0.0; 4 * n * 2 * n],
        }
    }

    #[inline]
    fn lead_to_idx(lead: i32) -> usize {
        if lead == 1 { 0 } else { 1 }
    }

    #[inline]
    fn rate_idx(n: usize, n1: usize, n2: usize, q1: usize, lead_idx: usize, q2: usize) -> usize {
        (((n1 * 2 + n2) * n + q1) * 2 + lead_idx) * n + q2
    }

    #[inline]
    pub fn rate_w(&self, n1: usize, n2: usize, q1: usize, q2: usize, lead: i32) -> f64 {
        if n1 > 1 || n2 > 1 || q1 >= self.n || q2 >= self.n {
            return 0.0;
        }
        let lead_idx = Self::lead_to_idx(lead);
        let idx = Self::rate_idx(self.n, n1, n2, q1, lead_idx, q2);
        self.data[idx]
    }

    #[inline]
    pub fn rate_w_lead(&self, n1: usize, n2: usize, q1: usize, q2: usize) -> f64 {
        self.rate_w(n1, n2, q1, q2, 1) + self.rate_w(n1, n2, q1, q2, -1)
    }
}

#[inline(always)]
fn sanitize(x: f64) -> f64 {
    if x.is_finite() { x } else { 0.0 }
}

pub fn spin_degeneracy(n1: usize, n2: usize) -> f64 {
    if n1 == 0 {
        2.0
    } else if n2 == 0 {
        1.0
    } else {
        2.0
    }
}

/// Compute rates for all q2 given (n1, n2, q1, lead).
pub fn m_rate_w(
    n1: usize,
    n2: usize,
    q1: i32,
    q2_vec: &[i32],
    vmode: f64,
    alpha_l: f64,
    alpha_r: f64,
    lambda: f64,
    vsd: f64,
    t: f64,
    eta: f64,
    lead: i32,
    vg: f64,
    fc: &mut FCCache,
    out: &mut [f64],
) {
    let nq2 = q2_vec.len();
    let gamma_l = alpha_l * vmode;
    let gamma_r = alpha_r * vmode;
    let epsilond = 0.0 + vg;
    let mu_l = eta * vsd;
    let mu_r = -(1.0 - eta) * vsd;

    let s = spin_degeneracy(n1, n2);
    let gamma = if lead == 1 { gamma_l } else { gamma_r };
    let mu = if lead == 1 { mu_l } else { mu_r };

    if lead != 1 && lead != -1 {
        for v in out.iter_mut().take(nq2) {
            *v = 0.0;
        }
        return;
    }

    if n1 == 1 && n2 == 0 {
        // Sequential 1→0
        for i in 0..nq2 {
            let q2 = q2_vec[i];
            let fc_val = fc.get(q1, q2);
            let f = fermi(epsilond - (q2 - q1) as f64 * vmode, mu, t);
            let val = s * gamma / HBAR_EV * fc_val * fc_val * (1.0 - f);
            out[i] = sanitize(val);
        }
    } else if n1 == 0 && n2 == 1 {
        // Sequential 0→1
        for i in 0..nq2 {
            let q2 = q2_vec[i];
            let fc_val = fc.get(q1, q2);
            let f = fermi(epsilond + (q2 - q1) as f64 * vmode, mu, t);
            let val = s * gamma / HBAR_EV * fc_val * fc_val * f;
            out[i] = sanitize(val);
        }
    } else if n1 == 0 && n2 == 0 {
        // Cotunneling 0→0
        let mut single_sum = vec![0.0; nq2];
        let mut double_sum = vec![0.0; nq2];

        if lead == 1 {
            sum_mmr(q1, q2_vec, lambda, mu_l, mu_r, vmode, epsilond, t, fc, &mut single_sum);
            sum_mmmmrs(q1, q2_vec, lambda, mu_l, mu_r, vmode, epsilond, t, fc, &mut double_sum);
        } else {
            sum_mmr(q1, q2_vec, lambda, mu_r, mu_l, vmode, epsilond, t, fc, &mut single_sum);
            sum_mmmmrs(q1, q2_vec, lambda, mu_r, mu_l, vmode, epsilond, t, fc, &mut double_sum);
        }

        let prefac = s / (2.0 * PI * HBAR_EV) * gamma_l * gamma_r;
        for i in 0..nq2 {
            out[i] = sanitize(prefac * (single_sum[i] + double_sum[i]));
        }
    } else {
        // Cotunneling 1→1
        let mut single_sum = vec![0.0; nq2];
        let mut double_sum = vec![0.0; nq2];

        if lead == 1 {
            sum_mmr11(q1, q2_vec, lambda, mu_l, mu_r, vmode, epsilond, t, fc, &mut single_sum);
            sum_mmmmrs11(
                q1, q2_vec, lambda, mu_l, mu_r, vmode, epsilond, t, fc, &mut double_sum,
            );
        } else {
            sum_mmr11(q1, q2_vec, lambda, mu_r, mu_l, vmode, epsilond, t, fc, &mut single_sum);
            sum_mmmmrs11(
                q1, q2_vec, lambda, mu_r, mu_l, vmode, epsilond, t, fc, &mut double_sum,
            );
        }

        let prefac = s / (2.0 * PI * HBAR_EV) * gamma_l * gamma_r;
        for i in 0..nq2 {
            out[i] = sanitize(prefac * (single_sum[i] + double_sum[i]));
        }
    }
}

/// Pre-compute all rates for a given lead and store them.
pub fn calculate_all_rate_w(
    store: &mut RateStore,
    n: usize,
    vmode: f64,
    alpha_l: f64,
    alpha_r: f64,
    lambda: f64,
    vsd: f64,
    t: f64,
    eta: f64,
    lead: i32,
    vg: f64,
    fc: &mut FCCache,
) {
    let q2_vec: Vec<i32> = (0..n as i32).collect();
    let mut out_buf = vec![0.0; n];
    let lead_idx = RateStore::lead_to_idx(lead);

    for n1 in 0..=1_usize {
        for n2 in 0..=1_usize {
            for q1 in 0..n {
                m_rate_w(
                    n1, n2, q1 as i32, &q2_vec, vmode, alpha_l, alpha_r, lambda, vsd, t, eta,
                    lead, vg, fc, &mut out_buf,
                );
                for q2 in 0..n {
                    let idx = RateStore::rate_idx(n, n1, n2, q1, lead_idx, q2);
                    store.data[idx] = out_buf[q2];
                }
            }
        }
    }
}
