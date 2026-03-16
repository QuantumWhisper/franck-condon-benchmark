use num_complex::Complex64;

use crate::digamma::{digamma_c, trigamma_c};

const TABLE_Y_MAX: f64 = 12.0;
const TABLE_DY: f64 = 0.0005;
const TABLE_SIZE: usize = (TABLE_Y_MAX / TABLE_DY) as usize + 1;
const DIRECT_THRESHOLD: f64 = 0.3;

pub struct DigammaTable {
    psi_re: Vec<f64>,
    psi_im: Vec<f64>,
    psi1_re: Vec<f64>,
    psi1_im: Vec<f64>,
}

impl DigammaTable {
    pub fn new() -> Self {
        let mut psi_re = vec![0.0; TABLE_SIZE];
        let mut psi_im = vec![0.0; TABLE_SIZE];
        let mut psi1_re = vec![0.0; TABLE_SIZE];
        let mut psi1_im = vec![0.0; TABLE_SIZE];

        for i in 0..TABLE_SIZE {
            let y = i as f64 * TABLE_DY;
            let z = Complex64::new(0.5, y);
            let psi = digamma_c(z);
            let psi1 = trigamma_c(z);
            psi_re[i] = psi.re;
            psi_im[i] = psi.im;
            psi1_re[i] = psi1.re;
            psi1_im[i] = psi1.im;
        }

        DigammaTable {
            psi_re,
            psi_im,
            psi1_re,
            psi1_im,
        }
    }

    #[inline(always)]
    fn catmull_rom(t: f64, p0: f64, p1: f64, p2: f64, p3: f64) -> f64 {
        let t2 = t * t;
        let t3 = t2 * t;
        0.5 * ((2.0 * p1)
            + (-p0 + p2) * t
            + (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t2
            + (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t3)
    }

    #[inline(always)]
    fn lookup_positive(&self, y: f64, re: &[f64], im: &[f64], use_trigamma_fallback: bool) -> Complex64 {
        let idx_f = y / TABLE_DY;
        let idx = idx_f as usize;

        if idx + 2 >= TABLE_SIZE {
            let z = Complex64::new(0.5, y);
            return if use_trigamma_fallback { trigamma_c(z) } else { digamma_c(z) };
        }

        let t = idx_f - idx as f64;

        let i0 = if idx > 0 { idx - 1 } else { 0 };
        let i1 = idx;
        let i2 = idx + 1;
        let i3 = idx + 2;

        let r = Self::catmull_rom(t, re[i0], re[i1], re[i2], re[i3]);
        let i = Self::catmull_rom(t, im[i0], im[i1], im[i2], im[i3]);
        Complex64::new(r, i)
    }

    #[inline]
    pub fn digamma(&self, y: f64) -> Complex64 {
        let abs_y = y.abs();

        if abs_y < DIRECT_THRESHOLD || abs_y >= TABLE_Y_MAX {
            return digamma_c(Complex64::new(0.5, y));
        }

        let val = self.lookup_positive(abs_y, &self.psi_re, &self.psi_im, false);

        if y >= 0.0 {
            val
        } else {
            val.conj()
        }
    }

    #[inline]
    pub fn trigamma(&self, y: f64) -> Complex64 {
        let abs_y = y.abs();

        if abs_y < DIRECT_THRESHOLD || abs_y >= TABLE_Y_MAX {
            return trigamma_c(Complex64::new(0.5, y));
        }

        let val = self.lookup_positive(abs_y, &self.psi1_re, &self.psi1_im, true);

        if y >= 0.0 {
            val
        } else {
            val.conj()
        }
    }
}

unsafe impl Send for DigammaTable {}
unsafe impl Sync for DigammaTable {}
