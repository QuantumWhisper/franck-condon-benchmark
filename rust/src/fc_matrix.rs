use crate::laguerre::laguerre_l;

pub const FC_MAX_N: usize = 256;

pub struct FCCache {
    cache: Vec<f64>,
    valid: Vec<bool>,
    pub lambda: f64,
}

impl FCCache {
    pub fn new(lambda: f64) -> Self {
        FCCache {
            cache: vec![0.0; FC_MAX_N * FC_MAX_N],
            valid: vec![false; FC_MAX_N * FC_MAX_N],
            lambda,
        }
    }

    #[inline]
    fn idx(q1: usize, q2: usize) -> usize {
        q1 * FC_MAX_N + q2
    }

    #[inline]
    pub fn get(&mut self, q1: i32, q2: i32) -> f64 {
        if q1 < 0 || q2 < 0 || q1 as usize >= FC_MAX_N || q2 as usize >= FC_MAX_N {
            return 0.0;
        }
        let i = Self::idx(q1 as usize, q2 as usize);
        if !self.valid[i] {
            self.cache[i] = fc_matrix_single(q1, q2, self.lambda);
            self.valid[i] = true;
        }
        self.cache[i]
    }

    pub fn get_row(&mut self, q1: i32, q2_range: &[i32], out: &mut [f64]) {
        for (i, &q2) in q2_range.iter().enumerate() {
            out[i] = self.get(q1, q2);
        }
    }
}

/// Compute log(q!) for non-negative integer q.
fn log_factorial(n: i32) -> f64 {
    if n <= 1 {
        return 0.0;
    }
    (2..=n).map(|k| (k as f64).ln()).sum()
}

/// Franck-Condon matrix element <q2|D(lambda)|q1>.
pub fn fc_matrix_single(q1: i32, q2: i32, lambda: f64) -> f64 {
    let q = q1.min(q2);
    let cap_q = q1.max(q2);

    if q1 == 0 && q2 == 0 {
        return (-lambda * lambda / 2.0).exp();
    }

    let mut sign_factor = 1.0;
    if q2 < q1 && ((q1 - q2) % 2 != 0) {
        sign_factor = -1.0;
    }

    let l = laguerre_l(q, cap_q - q, lambda * lambda);

    if lambda == 0.0 {
        return if q1 == q2 { l } else { 0.0 };
    }

    let log_coeff = (cap_q - q) as f64 * lambda.ln()
        - lambda * lambda / 2.0
        + 0.5 * (log_factorial(q) - log_factorial(cap_q));

    let m = sign_factor * log_coeff.exp() * l;
    if m.is_nan() || m.is_infinite() {
        0.0
    } else {
        m
    }
}
