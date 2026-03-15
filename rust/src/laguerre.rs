/// Generalized Laguerre polynomial L_n^alpha(x) via three-term recurrence.
pub fn laguerre_l(n: i32, alpha: i32, x: f64) -> f64 {
    if n <= 0 {
        return 1.0;
    }

    let mut l_prev = 1.0;
    let mut l_curr = 1.0 + alpha as f64 - x;
    if n == 1 {
        return l_curr;
    }

    for k in 1..n {
        let kf = k as f64;
        let af = alpha as f64;
        let l_next = ((2.0 * kf + 1.0 + af - x) * l_curr - (kf + af) * l_prev) / (kf + 1.0);
        l_prev = l_curr;
        l_curr = l_next;
    }

    l_curr
}
