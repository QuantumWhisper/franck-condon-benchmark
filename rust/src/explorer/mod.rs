use std::io;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{mpsc, Arc};
use std::thread;
use std::time::{Duration, Instant};

use ratatui::{
    crossterm::event::{self as crossterm_event, Event, KeyEventKind},
    DefaultTerminal,
};

use franck_condon::digamma_table::DigammaTable;
use franck_condon::fc_matrix::{FCCache, FC_MAX_N};
use franck_condon::simulate::{simulate_iv_with_cache, SimulationResult};

pub(crate) mod derivatives;
pub(crate) mod event;
pub(crate) mod export;
pub(crate) mod render;
pub(crate) mod types;
pub(crate) mod widgets;
pub(crate) use types::*;

const TICK_RATE: Duration = Duration::from_millis(50);
const PARAM_PANEL_WIDTH: u16 = 36;

// ---------------------------------------------------------------------------
// App state
// ---------------------------------------------------------------------------

pub(crate) struct App {
    pub(crate) mode: AppMode,
    pub(crate) should_quit: bool,
    pub(crate) selected_param: usize,

    // simulation parameters
    pub(crate) n: usize,
    pub(crate) lambda: f64,
    pub(crate) t_kelvin: f64,
    pub(crate) vmode: f64,
    pub(crate) alpha_l: f64,
    pub(crate) alpha_r: f64,
    pub(crate) eta: f64,
    pub(crate) vg: f64,
    pub(crate) tau: f64,

    // I-V sweep
    pub(crate) vsd_min: f64,
    pub(crate) vsd_max: f64,
    pub(crate) vsd_step: f64,

    // stability grid
    pub(crate) stab_vsd_max: f64,
    pub(crate) n_vsd: usize,
    pub(crate) vg_min: f64,
    pub(crate) vg_max: f64,
    pub(crate) n_vg: usize,

    // results
    pub(crate) iv_data: Option<(SimulationResult, Duration)>,
    pub(crate) stability_grid: Vec<Vec<f64>>,
    pub(crate) stability_vg_vals: Vec<f64>,
    pub(crate) stability_vsd_vals: Vec<f64>,
    pub(crate) stability_progress: (usize, usize),
    pub(crate) stability_elapsed: Option<Duration>,

    // Temperature diagram sweep
    pub(crate) temp_vsd_max: f64,
    pub(crate) temp_n_vsd: usize,
    pub(crate) t_min: f64,
    pub(crate) t_max: f64,
    pub(crate) n_t: usize,

    // Temperature diagram results
    pub(crate) temperature_grid: Vec<Vec<f64>>,
    pub(crate) temperature_t_vals: Vec<f64>,
    pub(crate) temperature_vsd_vals: Vec<f64>,
    pub(crate) temperature_progress: (usize, usize),
    pub(crate) temperature_elapsed: Option<Duration>,

    // background computation
    pub(crate) compute_start: Instant,
    pub(crate) compute_running: bool,
    pub(crate) compute_rx: mpsc::Receiver<ComputeMsg>,
    pub(crate) _keep_tx: mpsc::Sender<ComputeMsg>,
    pub(crate) cancel_flag: Arc<AtomicBool>,

    pub(crate) display_mode: DisplayMode,
    pub(crate) status_msg: Option<(String, bool)>, // (message, is_error)
    pub(crate) contrast_min: Option<f64>,
    pub(crate) contrast_max: Option<f64>,
    pub(crate) gamma: f64,
    pub(crate) cursor_mode: bool,
    pub(crate) cursor_col: usize,
    pub(crate) cursor_row: usize,
    pub(crate) show_line_profile: bool,
    pub(crate) line_profile_horizontal: bool,
    pub(crate) iv_recompute_at: Option<Instant>,
}

impl App {
    fn new() -> Self {
        let (tx, rx) = mpsc::channel();
        Self {
            mode: AppMode::IVCurve,
            should_quit: false,
            selected_param: 0,
            n: 6,
            lambda: 5.0,
            t_kelvin: 4.2,
            vmode: 0.073,
            alpha_l: 0.02,
            alpha_r: 0.02,
            eta: 0.5,
            vg: 0.0,
            tau: f64::INFINITY,
            vsd_min: 0.0,
            vsd_max: 0.6,
            vsd_step: 0.003,
            stab_vsd_max: 0.6,
            n_vsd: 101,
            vg_min: -0.3,
            vg_max: 0.3,
            n_vg: 61,
            iv_data: None,
            stability_grid: Vec::new(),
            stability_vg_vals: Vec::new(),
            stability_vsd_vals: Vec::new(),
            stability_progress: (0, 0),
            stability_elapsed: None,
            temp_vsd_max: 0.6,
            temp_n_vsd: 101,
            t_min: 1.0,
            t_max: 50.0,
            n_t: 50,
            temperature_grid: Vec::new(),
            temperature_t_vals: Vec::new(),
            temperature_vsd_vals: Vec::new(),
            temperature_progress: (0, 0),
            temperature_elapsed: None,
            compute_start: Instant::now(),
            compute_running: false,
            compute_rx: rx,
            _keep_tx: tx,
            cancel_flag: Arc::new(AtomicBool::new(false)),
            display_mode: DisplayMode::Current,
            status_msg: None,
            contrast_min: None,
            contrast_max: None,
            gamma: 1.0,
            cursor_mode: false,
            cursor_col: 0,
            cursor_row: 0,
            show_line_profile: false,
            line_profile_horizontal: true,
            iv_recompute_at: None,
        }
    }

    fn visible_params(&self) -> Vec<ParamId> {
        match self.mode {
            AppMode::IVCurve => vec![
                ParamId::N,
                ParamId::Lambda,
                ParamId::T,
                ParamId::Vmode,
                ParamId::AlphaL,
                ParamId::AlphaR,
                ParamId::Eta,
                ParamId::Vg,
                ParamId::Tau,
                ParamId::VsdMin,
                ParamId::VsdMax,
                ParamId::VsdStep,
            ],
            AppMode::Stability => vec![
                ParamId::N,
                ParamId::Lambda,
                ParamId::T,
                ParamId::Vmode,
                ParamId::AlphaL,
                ParamId::AlphaR,
                ParamId::Eta,
                ParamId::Tau,
                ParamId::StabVsdMax,
                ParamId::NVsd,
                ParamId::VgMin,
                ParamId::VgMax,
                ParamId::NVg,
            ],
            AppMode::Temperature => vec![
                ParamId::N,
                ParamId::Lambda,
                ParamId::Vmode,
                ParamId::AlphaL,
                ParamId::AlphaR,
                ParamId::Eta,
                ParamId::Vg,
                ParamId::Tau,
                ParamId::TempVsdMax,
                ParamId::TempNVsd,
                ParamId::TMin,
                ParamId::TMax,
                ParamId::NT,
            ],
        }
    }

    fn param_display(&self, id: ParamId) -> (&'static str, String) {
        match id {
            ParamId::N => ("N", format!("{}", self.n)),
            ParamId::Lambda => ("lambda", format!("{:.3}", self.lambda)),
            ParamId::T => ("T (K)", format!("{:.1}", self.t_kelvin)),
            ParamId::Vmode => ("hw (eV)", format!("{:.4}", self.vmode)),
            ParamId::AlphaL => ("alpha_L", format!("{:.4}", self.alpha_l)),
            ParamId::AlphaR => ("alpha_R", format!("{:.4}", self.alpha_r)),
            ParamId::Eta => ("eta", format!("{:.2}", self.eta)),
            ParamId::Vg => ("Vg (V)", format!("{:.3}", self.vg)),
            ParamId::Tau => (
                "tau",
                if self.tau.is_infinite() {
                    "Inf".into()
                } else {
                    format!("{:.2e}", self.tau)
                },
            ),
            ParamId::VsdMin => ("Vsd min", format!("{:.3}", self.vsd_min)),
            ParamId::VsdMax => ("Vsd max", format!("{:.3}", self.vsd_max)),
            ParamId::VsdStep => ("Vsd step", format!("{:.4}", self.vsd_step)),
            ParamId::StabVsdMax => ("Vsd max", format!("{:.3}", self.stab_vsd_max)),
            ParamId::NVsd => ("# Vsd", format!("{}", self.n_vsd)),
            ParamId::VgMin => ("Vg min", format!("{:.3}", self.vg_min)),
            ParamId::VgMax => ("Vg max", format!("{:.3}", self.vg_max)),
            ParamId::NVg => ("# Vg", format!("{}", self.n_vg)),
            ParamId::TempVsdMax => ("Vsd max", format!("{:.3}", self.temp_vsd_max)),
            ParamId::TempNVsd => ("# Vsd", format!("{}", self.temp_n_vsd)),
            ParamId::TMin => ("T min (K)", format!("{:.1}", self.t_min)),
            ParamId::TMax => ("T max (K)", format!("{:.1}", self.t_max)),
            ParamId::NT => ("# T", format!("{}", self.n_t)),
        }
    }

    fn adjust_param(&mut self, id: ParamId, dir: i32, fine: bool) {
        let d = dir as f64;
        match id {
            ParamId::N => {
                self.n = (self.n as i32 + dir).clamp(1, 200) as usize;
            }
            ParamId::Lambda => {
                let s = if fine { 0.1 } else { 0.5 };
                self.lambda = (self.lambda + d * s).clamp(0.0, 15.0);
            }
            ParamId::T => {
                let s = if fine { 0.1 } else { 1.0 };
                self.t_kelvin = (self.t_kelvin + d * s).max(0.1);
            }
            ParamId::Vmode => {
                let s = if fine { 0.001 } else { 0.01 };
                self.vmode = (self.vmode + d * s).max(0.001);
            }
            ParamId::AlphaL => {
                let s = if fine { 0.001 } else { 0.01 };
                self.alpha_l = (self.alpha_l + d * s).max(0.0);
            }
            ParamId::AlphaR => {
                let s = if fine { 0.001 } else { 0.01 };
                self.alpha_r = (self.alpha_r + d * s).max(0.0);
            }
            ParamId::Eta => {
                let s = if fine { 0.01 } else { 0.1 };
                self.eta = (self.eta + d * s).clamp(0.0, 1.0);
            }
            ParamId::Vg => {
                let s = if fine { 0.01 } else { 0.1 };
                self.vg += d * s;
            }
            ParamId::Tau => {
                if self.tau.is_infinite() {
                    if dir < 0 {
                        self.tau = 1e-9;
                    }
                } else if dir > 0 {
                    self.tau *= 10.0;
                    if self.tau > 1e6 {
                        self.tau = f64::INFINITY;
                    }
                } else {
                    self.tau /= 10.0;
                    if self.tau < 1e-15 {
                        self.tau = 1e-15;
                    }
                }
            }
            ParamId::VsdMin => {
                let s = if fine { 0.01 } else { 0.1 };
                self.vsd_min = (self.vsd_min + d * s).max(0.0);
            }
            ParamId::VsdMax => {
                let s = if fine { 0.01 } else { 0.1 };
                self.vsd_max = (self.vsd_max + d * s).max(0.01);
            }
            ParamId::VsdStep => {
                let s = if fine { 0.0005 } else { 0.001 };
                self.vsd_step = (self.vsd_step + d * s).clamp(0.001, 0.1);
            }
            ParamId::StabVsdMax => {
                let s = if fine { 0.01 } else { 0.1 };
                self.stab_vsd_max = (self.stab_vsd_max + d * s).max(0.05);
            }
            ParamId::NVsd => {
                let step = if fine { 1 } else { 10 };
                self.n_vsd = (self.n_vsd as i32 + dir * step).clamp(10, 501) as usize;
            }
            ParamId::VgMin => {
                let s = if fine { 0.01 } else { 0.1 };
                self.vg_min += d * s;
            }
            ParamId::VgMax => {
                let s = if fine { 0.01 } else { 0.1 };
                self.vg_max += d * s;
            }
            ParamId::NVg => {
                let step = if fine { 1 } else { 10 };
                self.n_vg = (self.n_vg as i32 + dir * step).clamp(10, 501) as usize;
            }
            ParamId::TempVsdMax => {
                let s = if fine { 0.01 } else { 0.1 };
                self.temp_vsd_max = (self.temp_vsd_max + d * s).max(0.05);
            }
            ParamId::TempNVsd => {
                let step = if fine { 1 } else { 10 };
                self.temp_n_vsd = (self.temp_n_vsd as i32 + dir * step).clamp(10, 501) as usize;
            }
            ParamId::TMin => {
                let s = if fine { 0.1 } else { 1.0 };
                self.t_min = (self.t_min + d * s).max(0.1);
            }
            ParamId::TMax => {
                let s = if fine { 0.1 } else { 5.0 };
                self.t_max = (self.t_max + d * s).max(0.5);
            }
            ParamId::NT => {
                let step = if fine { 1 } else { 10 };
                self.n_t = (self.n_t as i32 + dir * step).clamp(5, 501) as usize;
            }
        }
    }

    fn validate(&self) -> Result<(), String> {
        if self.t_kelvin <= 0.0 {
            return Err("T must be > 0".into());
        }
        if self.vmode <= 0.0 {
            return Err("hw must be > 0".into());
        }
        if self.n == 0 {
            return Err("N must be >= 1".into());
        }
        if self.lambda < 0.0 {
            return Err("lambda must be >= 0".into());
        }
        if self.alpha_l < 0.0 || self.alpha_r < 0.0 {
            return Err("alpha must be >= 0".into());
        }
        if !self.tau.is_infinite() && self.tau <= 0.0 {
            return Err("tau must be > 0 or Inf".into());
        }
        Ok(())
    }

    fn fc_precompute_bound(&self) -> usize {
        ((self.lambda.powf(2.2) * 3.0).round() as usize + 50)
            .max(self.n)
            .min(FC_MAX_N)
    }

    fn start_iv(&mut self) {
        if self.compute_running {
            return;
        }
        if let Err(e) = self.validate() {
            self.status_msg = Some((e, true));
            return;
        }
        self.status_msg = None;

        let (tx, rx) = mpsc::channel();
        self.compute_rx = rx;
        self.cancel_flag = Arc::new(AtomicBool::new(false));
        self.compute_start = Instant::now();
        self.compute_running = true;

        let n = self.n;
        let vmode = self.vmode;
        let al = self.alpha_l;
        let ar = self.alpha_r;
        let lambda = self.lambda;
        let t = self.t_kelvin;
        let eta = self.eta;
        let vg = self.vg;
        let tau = self.tau;
        let vsd_min = self.vsd_min;
        let vsd_max = self.vsd_max;
        let vsd_step = self.vsd_step;
        let bound = self.fc_precompute_bound();

        thread::spawn(move || {
            let start = Instant::now();

            let n_pts = ((vsd_max - vsd_min) / vsd_step).round() as usize + 1;
            let vsd_vec: Vec<f64> = (0..n_pts).map(|i| vsd_min + i as f64 * vsd_step).collect();

            let mut fc = FCCache::new(lambda);
            fc.precompute(bound);
            let dtable = DigammaTable::new();

            let result = simulate_iv_with_cache(
                n, vmode, al, ar, lambda, &vsd_vec, t, eta, vg, tau, &fc, &dtable,
            );

            let elapsed = start.elapsed();
            tx.send(ComputeMsg::IVDone { result, elapsed }).ok();
        });
    }

    fn start_stability(&mut self) {
        if self.compute_running {
            return;
        }
        if let Err(e) = self.validate() {
            self.status_msg = Some((e, true));
            return;
        }
        self.status_msg = None;

        let (tx, rx) = mpsc::channel();
        self.compute_rx = rx;
        self.cancel_flag = Arc::new(AtomicBool::new(false));
        self.compute_start = Instant::now();
        self.compute_running = true;

        let vg_vals: Vec<f64> = (0..self.n_vg)
            .map(|i| {
                self.vg_min + i as f64 * (self.vg_max - self.vg_min) / (self.n_vg - 1).max(1) as f64
            })
            .collect();
        let vsd_vals: Vec<f64> = (0..self.n_vsd)
            .map(|i| i as f64 * self.stab_vsd_max / (self.n_vsd - 1).max(1) as f64)
            .collect();

        self.stability_grid = vec![Vec::new(); vg_vals.len()];
        self.stability_vg_vals = vg_vals.clone();
        self.stability_vsd_vals = vsd_vals.clone();
        self.stability_progress = (0, vg_vals.len());
        self.stability_elapsed = None;

        let n = self.n;
        let vmode = self.vmode;
        let al = self.alpha_l;
        let ar = self.alpha_r;
        let lambda = self.lambda;
        let t = self.t_kelvin;
        let eta = self.eta;
        let tau = self.tau;
        let cancel = self.cancel_flag.clone();
        let bound = self.fc_precompute_bound();

        thread::spawn(move || {
            let start = Instant::now();

            let mut fc = FCCache::new(lambda);
            fc.precompute(bound);
            let dtable = DigammaTable::new();

            for (i, &vg) in vg_vals.iter().enumerate() {
                if cancel.load(Ordering::Relaxed) {
                    break;
                }
                let result = simulate_iv_with_cache(
                    n, vmode, al, ar, lambda, &vsd_vals, t, eta, vg, tau, &fc, &dtable,
                );
                if tx
                    .send(ComputeMsg::StabilityRow {
                        vg_idx: i,
                        i_tol: result.i_tol,
                    })
                    .is_err()
                {
                    break;
                }
            }
            tx.send(ComputeMsg::StabilityDone {
                elapsed: start.elapsed(),
            })
            .ok();
        });
    }

    fn start_temperature(&mut self) {
        if self.compute_running {
            return;
        }
        if self.t_min <= 0.0 {
            self.status_msg = Some(("T min must be > 0".into(), true));
            return;
        }
        if self.t_min >= self.t_max {
            self.status_msg = Some(("T min must be < T max".into(), true));
            return;
        }
        if let Err(e) = self.validate() {
            self.status_msg = Some((e, true));
            return;
        }
        self.status_msg = None;

        let (tx, rx) = mpsc::channel();
        self.compute_rx = rx;
        self.cancel_flag = Arc::new(AtomicBool::new(false));
        self.compute_start = Instant::now();
        self.compute_running = true;

        let t_vals: Vec<f64> = (0..self.n_t)
            .map(|i| {
                self.t_min + i as f64 * (self.t_max - self.t_min) / (self.n_t - 1).max(1) as f64
            })
            .collect();
        let vsd_vals: Vec<f64> = (0..self.temp_n_vsd)
            .map(|i| i as f64 * self.temp_vsd_max / (self.temp_n_vsd - 1).max(1) as f64)
            .collect();

        self.temperature_grid = vec![Vec::new(); t_vals.len()];
        self.temperature_t_vals = t_vals.clone();
        self.temperature_vsd_vals = vsd_vals.clone();
        self.temperature_progress = (0, t_vals.len());
        self.temperature_elapsed = None;

        let n = self.n;
        let vmode = self.vmode;
        let al = self.alpha_l;
        let ar = self.alpha_r;
        let lambda = self.lambda;
        let eta = self.eta;
        let vg = self.vg;
        let tau = self.tau;
        let cancel = self.cancel_flag.clone();
        let bound = self.fc_precompute_bound();

        thread::spawn(move || {
            let start = Instant::now();

            let mut fc = FCCache::new(lambda);
            fc.precompute(bound);
            let dtable = DigammaTable::new();

            for (i, &t) in t_vals.iter().enumerate() {
                if cancel.load(Ordering::Relaxed) {
                    break;
                }
                let result = simulate_iv_with_cache(
                    n, vmode, al, ar, lambda, &vsd_vals, t, eta, vg, tau, &fc, &dtable,
                );
                if tx
                    .send(ComputeMsg::TemperatureRow {
                        t_idx: i,
                        i_tol: result.i_tol,
                    })
                    .is_err()
                {
                    break;
                }
            }
            tx.send(ComputeMsg::TemperatureDone {
                elapsed: start.elapsed(),
            })
            .ok();
        });
    }

    fn cancel_compute(&mut self) {
        self.cancel_flag.store(true, Ordering::Relaxed);
        self.compute_running = false;
    }

    fn poll_messages(&mut self) {
        while let Ok(msg) = self.compute_rx.try_recv() {
            match msg {
                ComputeMsg::IVDone { result, elapsed } => {
                    self.status_msg =
                        Some((format!("Done ({:.2}s)", elapsed.as_secs_f64()), false));
                    self.iv_data = Some((result, elapsed));
                    self.compute_running = false;
                }
                ComputeMsg::StabilityRow { vg_idx, i_tol } => {
                    if vg_idx < self.stability_grid.len() {
                        self.stability_grid[vg_idx] = i_tol;
                        self.stability_progress.0 = vg_idx + 1;
                    }
                }
                ComputeMsg::StabilityDone { elapsed } => {
                    self.status_msg =
                        Some((format!("Done ({:.1}s)", elapsed.as_secs_f64()), false));
                    self.stability_elapsed = Some(elapsed);
                    self.compute_running = false;
                }
                ComputeMsg::TemperatureRow { t_idx, i_tol } => {
                    if t_idx < self.temperature_grid.len() {
                        self.temperature_grid[t_idx] = i_tol;
                        self.temperature_progress.0 = t_idx + 1;
                    }
                }
                ComputeMsg::TemperatureDone { elapsed } => {
                    self.status_msg =
                        Some((format!("Done ({:.1}s)", elapsed.as_secs_f64()), false));
                    self.temperature_elapsed = Some(elapsed);
                    self.compute_running = false;
                }
            }
        }

        if let Some(deadline) = self.iv_recompute_at {
            if Instant::now() >= deadline && !self.compute_running {
                self.iv_recompute_at = None;
                self.start_iv();
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

fn main() -> io::Result<()> {
    let terminal = ratatui::init();
    let result = run(terminal);
    ratatui::restore();
    result
}

fn run(mut terminal: DefaultTerminal) -> io::Result<()> {
    let mut app = App::new();
    let mut last_tick = Instant::now();

    loop {
        terminal.draw(|frame| render::render(frame, &app))?;

        let timeout = TICK_RATE.saturating_sub(last_tick.elapsed());
        if crossterm_event::poll(timeout)? {
            if let Event::Key(key) = crossterm_event::read()? {
                if key.kind == KeyEventKind::Press {
                    event::handle_key(&mut app, key);
                }
            }
        }

        if last_tick.elapsed() >= TICK_RATE {
            app.poll_messages();
            last_tick = Instant::now();
        }

        if app.should_quit {
            break;
        }
    }

    Ok(())
}
