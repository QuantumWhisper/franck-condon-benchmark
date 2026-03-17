pub(crate) fn compute_didv(vsd: &[f64], current: &[f64]) -> Vec<f64> {
    let n = vsd.len().min(current.len());
    if n < 2 {
        return vec![0.0; n];
    }
    let mut g = vec![0.0; n];
    let dv0 = vsd[1] - vsd[0];
    g[0] = if dv0.abs() > 1e-30 {
        (current[1] - current[0]) / dv0
    } else {
        0.0
    };
    for i in 1..n - 1 {
        let dv = vsd[i + 1] - vsd[i - 1];
        g[i] = if dv.abs() > 1e-30 {
            (current[i + 1] - current[i - 1]) / dv
        } else {
            0.0
        };
    }
    let dvn = vsd[n - 1] - vsd[n - 2];
    g[n - 1] = if dvn.abs() > 1e-30 {
        (current[n - 1] - current[n - 2]) / dvn
    } else {
        0.0
    };
    g
}

pub(crate) fn compute_didv_grid(grid: &[Vec<f64>], vsd_vals: &[f64]) -> Vec<Vec<f64>> {
    grid.iter()
        .map(|row| {
            if row.is_empty() {
                Vec::new()
            } else {
                compute_didv(vsd_vals, row)
            }
        })
        .collect()
}

pub(crate) fn compute_d2idv2(vsd: &[f64], current: &[f64]) -> Vec<f64> {
    let n = vsd.len().min(current.len());
    if n < 3 {
        return vec![0.0; n];
    }
    let mut d2 = vec![0.0; n];
    for i in 1..n - 1 {
        let h1 = vsd[i] - vsd[i - 1];
        let h2 = vsd[i + 1] - vsd[i];
        let denom = h1 * h2 * (h1 + h2) / 2.0;
        d2[i] = if denom.abs() > 1e-60 {
            (current[i - 1] * h2 - current[i] * (h1 + h2) + current[i + 1] * h1) / denom
        } else {
            0.0
        };
    }
    d2[0] = d2[1];
    d2[n - 1] = d2[n - 2];
    d2
}

pub(crate) fn compute_d2idv2_grid(grid: &[Vec<f64>], vsd_vals: &[f64]) -> Vec<Vec<f64>> {
    grid.iter()
        .map(|row| {
            if row.is_empty() {
                Vec::new()
            } else {
                compute_d2idv2(vsd_vals, row)
            }
        })
        .collect()
}

pub(crate) fn compute_normalized_iets(vsd: &[f64], current: &[f64]) -> Vec<f64> {
    let didv = compute_didv(vsd, current);
    let d2idv2 = compute_d2idv2(vsd, current);
    didv.iter()
        .zip(d2idv2.iter())
        .map(|(&g, &d2)| if g.abs() > 1e-30 { d2 / g } else { 0.0 })
        .collect()
}

pub(crate) fn compute_normalized_iets_grid(grid: &[Vec<f64>], vsd_vals: &[f64]) -> Vec<Vec<f64>> {
    grid.iter()
        .map(|row| {
            if row.is_empty() {
                Vec::new()
            } else {
                compute_normalized_iets(vsd_vals, row)
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_didv_linear() {
        let vsd = vec![0.0, 1.0, 2.0, 3.0];
        let current = vec![0.0, 2.0, 4.0, 6.0];
        let g = compute_didv(&vsd, &current);
        for &val in &g {
            assert!((val - 2.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_d2idv2_quadratic() {
        let vsd = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let current = vec![0.0, 1.0, 4.0, 9.0, 16.0];
        let d2 = compute_d2idv2(&vsd, &current);
        for val in d2.iter().take(4).skip(1) {
            assert!((*val - 2.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_normalized_iets_zero_guard() {
        let vsd = vec![0.0, 1.0, 2.0];
        let current = vec![0.0, 0.0, 0.0];
        let niets = compute_normalized_iets(&vsd, &current);
        for &val in &niets {
            assert!(val.is_finite());
        }
    }

    #[test]
    fn test_didv_empty() {
        assert!(compute_didv(&[], &[]).is_empty());
    }

    #[test]
    fn test_didv_single() {
        let r = compute_didv(&[1.0], &[5.0]);
        assert_eq!(r.len(), 1);
        assert_eq!(r[0], 0.0);
    }

    #[test]
    fn test_d2idv2_short() {
        let r = compute_d2idv2(&[0.0, 1.0], &[0.0, 1.0]);
        assert_eq!(r.len(), 2);
    }

    #[test]
    fn test_grid_with_empty_rows() {
        let grid = vec![vec![1.0, 2.0], vec![], vec![3.0, 4.0]];
        let vsd = vec![0.0, 1.0];
        let result = compute_didv_grid(&grid, &vsd);
        assert_eq!(result.len(), 3);
        assert!(result[1].is_empty());
    }
}
