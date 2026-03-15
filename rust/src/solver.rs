use nalgebra::{DMatrix, DVector};

/// Solve steady-state occupation probabilities from rate equation 0 = WP.
/// Uses augmented system: replace last row with normalization sum(P)=1,
/// solve via QR decomposition, clamp negatives, renormalize.
pub fn solve_steady_state(w: &[f64], n: usize) -> Vec<f64> {
    let dim = 2 * n;
    if dim == 0 {
        return vec![];
    }

    // Build augmented matrix: copy W, replace last row with 1s
    let a = DMatrix::from_fn(dim, dim, |i, j| {
        if i == dim - 1 {
            1.0
        } else {
            w[i * dim + j]
        }
    });

    let mut d = DVector::zeros(dim);
    d[dim - 1] = 1.0;

    // Solve via QR decomposition
    let qr = a.qr();
    let x = match qr.solve(&d) {
        Some(sol) => sol,
        None => return vec![0.0; dim],
    };

    // Clamp negatives, renormalize
    let mut p: Vec<f64> = x
        .iter()
        .map(|&v| if v.is_finite() && v >= 0.0 { v } else { 0.0 })
        .collect();

    let sum: f64 = p.iter().sum();
    if sum > 0.0 {
        for v in &mut p {
            *v /= sum;
        }
    } else {
        p = vec![0.0; dim];
        p[dim - 1] = 1.0;
    }

    p
}
