use crate::constants::KB_EV;
use crate::rate::RateStore;

/// Equilibrium phonon occupation probability.
#[inline]
pub fn peq(q: i32, vmode: f64, t: f64) -> f64 {
    let beta = 1.0 / (KB_EV * t);
    (-(q as f64) * vmode * beta).exp() * (1.0 - (-vmode * beta).exp())
}

/// Sum of rates over all final phonon states.
fn sigma_w(store: &RateStore, n: usize, n1: usize, n2: usize, q1: usize) -> f64 {
    let mut s = 0.0;
    for q2 in 0..n {
        s += store.rate_w_lead(n1, n2, q1, q2);
    }
    s
}

/// Generate the 2N x 2N rate equation matrix W.
/// Output is stored in row-major order: w[i * dim + j].
pub fn generate_matrix_w(
    store: &RateStore,
    n: usize,
    vmode: f64,
    t: f64,
    tau: f64,
    w: &mut [f64],
) {
    let dim = 2 * n;
    let inv_tau = 1.0 / tau; // 1/INFINITY == 0 per IEEE 754

    // First block: n=0 rows (rows 0..N-1)
    for ii in 0..n {
        for jj in 0..dim {
            let modjj = jj % n;
            let p_q_index = modjj;
            let val;

            if jj < n {
                if ii == modjj {
                    // Diagonal of n=0 block
                    val = store.rate_w_lead(0, 0, p_q_index, p_q_index)
                        - sigma_w(store, n, 0, 0, p_q_index)
                        - sigma_w(store, n, 0, 1, p_q_index)
                        - inv_tau
                        + peq(p_q_index as i32, vmode, t) * inv_tau;
                } else {
                    // Off-diagonal of n=0 block
                    val = store.rate_w_lead(0, 0, p_q_index, ii)
                        + peq(ii as i32, vmode, t) * inv_tau;
                }
            } else {
                // n=1 contribution to n=0 rows
                val = store.rate_w_lead(1, 0, p_q_index, ii);
            }

            w[ii * dim + jj] = val;
        }
    }

    // Second block: n=1 rows (rows N..2N-1)
    for ii in 0..n {
        let row = n + ii;
        for jj in 0..dim {
            let modjj = jj % n;
            let val;

            if jj < n {
                val = store.rate_w_lead(0, 1, modjj, ii);
            } else {
                if modjj == ii {
                    // Diagonal of n=1 block
                    val = store.rate_w_lead(1, 1, modjj, ii)
                        - sigma_w(store, n, 1, 0, ii)
                        - sigma_w(store, n, 1, 1, ii)
                        - inv_tau
                        + peq(ii as i32, vmode, t) * inv_tau;
                } else {
                    // Off-diagonal of n=1 block
                    val = store.rate_w_lead(1, 1, modjj, ii)
                        + peq(ii as i32, vmode, t) * inv_tau;
                }
            }

            w[row * dim + jj] = val;
        }
    }
}
