use std::io::{self, Write};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{mpsc, Arc};
use std::thread;
use std::time::{Duration, Instant};

use ratatui::{
    buffer::Buffer,
    crossterm::event::{self, Event, KeyCode, KeyEventKind, KeyModifiers},
    layout::{Alignment, Constraint, Layout, Rect},
    style::{Color, Modifier, Style},
    symbols,
    text::Line,
    widgets::{Axis, Block, Chart, Dataset, Gauge, GraphType, LegendPosition, Paragraph, Widget},
    DefaultTerminal, Frame,
};

use franck_condon::digamma_table::DigammaTable;
use franck_condon::fc_matrix::{FCCache, FC_MAX_N};
use franck_condon::simulate::{simulate_iv_with_cache, SimulationResult};

const TICK_RATE: Duration = Duration::from_millis(50);
const PARAM_PANEL_WIDTH: u16 = 36;

// Viridis colormap — 9 control points from the official matplotlib colormap
const VIRIDIS: [(f64, f64, f64); 9] = [
    (0.267, 0.004, 0.329),
    (0.282, 0.140, 0.458),
    (0.253, 0.265, 0.530),
    (0.206, 0.372, 0.553),
    (0.127, 0.566, 0.551),
    (0.206, 0.718, 0.473),
    (0.430, 0.801, 0.348),
    (0.707, 0.868, 0.173),
    (0.993, 0.906, 0.144),
];

fn viridis(t: f64) -> Color {
    let t = t.clamp(0.0, 1.0);
    let n = VIRIDIS.len() - 1;
    let scaled = t * n as f64;
    let i = (scaled as usize).min(n - 1);
    let f = scaled - i as f64;
    let (r0, g0, b0) = VIRIDIS[i];
    let (r1, g1, b1) = VIRIDIS[i + 1];
    Color::Rgb(
        ((r0 + (r1 - r0) * f) * 255.0) as u8,
        ((g0 + (g1 - g0) * f) * 255.0) as u8,
        ((b0 + (b1 - b0) * f) * 255.0) as u8,
    )
}

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq, Eq)]
enum AppMode {
    IVCurve,
    Stability,
}

impl AppMode {
    fn label(self) -> &'static str {
        match self {
            Self::IVCurve => "I-V Curve",
            Self::Stability => "Stability Diagram",
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum ParamId {
    N,
    Lambda,
    T,
    Vmode,
    AlphaL,
    AlphaR,
    Eta,
    Vg,
    Tau,
    VsdMin,
    VsdMax,
    VsdStep,
    StabVsdMax,
    NVsd,
    VgMin,
    VgMax,
    NVg,
}

enum ComputeMsg {
    IVDone {
        result: SimulationResult,
        elapsed: Duration,
    },
    StabilityRow {
        vg_idx: usize,
        i_tol: Vec<f64>,
    },
    StabilityDone {
        elapsed: Duration,
    },
}

// ---------------------------------------------------------------------------
// App state
// ---------------------------------------------------------------------------

struct App {
    mode: AppMode,
    should_quit: bool,
    selected_param: usize,

    // simulation parameters
    n: usize,
    lambda: f64,
    t_kelvin: f64,
    vmode: f64,
    alpha_l: f64,
    alpha_r: f64,
    eta: f64,
    vg: f64,
    tau: f64,

    // I-V sweep
    vsd_min: f64,
    vsd_max: f64,
    vsd_step: f64,

    // stability grid
    stab_vsd_max: f64,
    n_vsd: usize,
    vg_min: f64,
    vg_max: f64,
    n_vg: usize,

    // results
    iv_data: Option<(SimulationResult, Duration)>,
    stability_grid: Vec<Vec<f64>>,
    stability_vg_vals: Vec<f64>,
    stability_vsd_vals: Vec<f64>,
    stability_progress: (usize, usize),
    stability_elapsed: Option<Duration>,

    // background computation
    compute_start: Instant,
    compute_running: bool,
    compute_rx: mpsc::Receiver<ComputeMsg>,
    _keep_tx: mpsc::Sender<ComputeMsg>,
    cancel_flag: Arc<AtomicBool>,

    status_msg: Option<(String, bool)>, // (message, is_error)
    iv_recompute_at: Option<Instant>,
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
            compute_start: Instant::now(),
            compute_running: false,
            compute_rx: rx,
            _keep_tx: tx,
            cancel_flag: Arc::new(AtomicBool::new(false)),
            status_msg: None,
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
            let vsd_vec: Vec<f64> = (0..n_pts)
                .map(|i| vsd_min + i as f64 * vsd_step)
                .collect();

            let mut fc = FCCache::new(lambda);
            fc.precompute(bound);
            let dtable = DigammaTable::new();

            let result =
                simulate_iv_with_cache(n, vmode, al, ar, lambda, &vsd_vec, t, eta, vg, tau, &fc, &dtable);

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
                self.vg_min
                    + i as f64 * (self.vg_max - self.vg_min) / (self.n_vg - 1).max(1) as f64
            })
            .collect();
        let vsd_vals: Vec<f64> = (0..self.n_vsd)
            .map(|i| {
                i as f64 * self.stab_vsd_max / (self.n_vsd - 1).max(1) as f64
            })
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

    fn cancel_compute(&mut self) {
        self.cancel_flag.store(true, Ordering::Relaxed);
        self.compute_running = false;
    }

    fn poll_messages(&mut self) {
        while let Ok(msg) = self.compute_rx.try_recv() {
            match msg {
                ComputeMsg::IVDone { result, elapsed } => {
                    self.status_msg = Some((format!("Done ({:.2}s)", elapsed.as_secs_f64()), false));
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
            }
        }

        if let Some(deadline) = self.iv_recompute_at {
            if Instant::now() >= deadline && !self.compute_running {
                self.iv_recompute_at = None;
                self.start_iv();
            }
        }
    }

    fn export_stability(&self) -> io::Result<String> {
        let filename = "stability_export.csv";
        let mut f = std::fs::File::create(filename)?;
        writeln!(f, "Vg_V,Vsd_V,I_tol_A")?;
        for (i, &vg) in self.stability_vg_vals.iter().enumerate() {
            if i >= self.stability_grid.len() || self.stability_grid[i].is_empty() {
                continue;
            }
            for (j, &vsd) in self.stability_vsd_vals.iter().enumerate() {
                if j < self.stability_grid[i].len() {
                    writeln!(f, "{:.6e},{:.6e},{:.6e}", vg, vsd, self.stability_grid[i][j])?;
                }
            }
        }
        Ok(filename.into())
    }

    fn export_iv(&self) -> io::Result<String> {
        let filename = "iv_export.csv";
        if let Some((ref data, _)) = self.iv_data {
            let mut f = std::fs::File::create(filename)?;
            writeln!(f, "Vsd_V,I_tol_A,I_seq_A,I_cot_A")?;
            for i in 0..data.vsd.len() {
                writeln!(
                    f,
                    "{:.6e},{:.6e},{:.6e},{:.6e}",
                    data.vsd[i], data.i_tol[i], data.i_seq[i], data.i_cot[i]
                )?;
            }
            Ok(filename.into())
        } else {
            Err(io::Error::new(io::ErrorKind::Other, "No I-V data"))
        }
    }
}

// ---------------------------------------------------------------------------
// UI rendering
// ---------------------------------------------------------------------------

fn render(frame: &mut Frame, app: &App) {
    let [main_area, help_area] =
        Layout::vertical([Constraint::Fill(1), Constraint::Length(1)]).areas(frame.area());

    let [params_area, plot_area] =
        Layout::horizontal([Constraint::Length(PARAM_PANEL_WIDTH), Constraint::Fill(1)])
            .areas(main_area);

    render_params(frame, params_area, app);

    match app.mode {
        AppMode::IVCurve => render_iv_chart(frame, plot_area, app),
        AppMode::Stability => render_stability(frame, plot_area, app),
    }

    render_help(frame, help_area, app);
}

fn render_params(frame: &mut Frame, area: Rect, app: &App) {
    let params = app.visible_params();
    let block = Block::bordered().title(Line::from(format!(" {} ", app.mode.label())).centered());
    let inner = block.inner(area);
    frame.render_widget(block, area);

    for (i, &id) in params.iter().enumerate() {
        if i as u16 >= inner.height {
            break;
        }
        let (name, value) = app.param_display(id);
        let selected = i == app.selected_param;

        let row = Rect {
            x: inner.x,
            y: inner.y + i as u16,
            width: inner.width,
            height: 1,
        };

        let prefix = if selected { " > " } else { "   " };
        let style = if selected {
            Style::new()
                .fg(Color::Yellow)
                .add_modifier(Modifier::BOLD)
        } else {
            Style::new().fg(Color::Gray)
        };

        let text = format!("{}{:<12}{:>10}", prefix, name, value);
        frame.render_widget(Paragraph::new(text).style(style), row);
    }

    let status_y = inner.y + (params.len() as u16).min(inner.height) + 1;
    if status_y < inner.bottom().saturating_sub(1) {
        let status_area = Rect {
            x: inner.x,
            y: status_y,
            width: inner.width,
            height: inner.bottom().saturating_sub(status_y),
        };

        let (text, style) = if app.compute_running {
            let dots = ".".repeat((app.compute_start.elapsed().as_millis() as usize / 500 % 4) + 1);
            (
                format!("  Computing{:<3}", dots),
                Style::new().fg(Color::Cyan),
            )
        } else if let Some((ref msg, is_err)) = app.status_msg {
            let c = if is_err { Color::Red } else { Color::Green };
            (format!("  {}", msg), Style::new().fg(c))
        } else {
            (
                "  Press Enter".into(),
                Style::new().fg(Color::DarkGray),
            )
        };

        frame.render_widget(Paragraph::new(text).style(style), status_area);
    }
}

fn render_iv_chart(frame: &mut Frame, area: Rect, app: &App) {
    let title = format!(
        " I-V Curve  N={}  l={:.1}  T={:.1}K ",
        app.n, app.lambda, app.t_kelvin
    );
    let block = Block::bordered().title(Line::from(title).centered());

    if let Some((ref data, _)) = app.iv_data {
        let data_tol: Vec<(f64, f64)> = data
            .vsd
            .iter()
            .zip(data.i_tol.iter())
            .map(|(&x, &y)| (x, y))
            .collect();
        let data_seq: Vec<(f64, f64)> = data
            .vsd
            .iter()
            .zip(data.i_seq.iter())
            .map(|(&x, &y)| (x, y))
            .collect();
        let data_cot: Vec<(f64, f64)> = data
            .vsd
            .iter()
            .zip(data.i_cot.iter())
            .map(|(&x, &y)| (x, y))
            .collect();

        let x_min = data.vsd.first().copied().unwrap_or(0.0);
        let x_max = data.vsd.last().copied().unwrap_or(1.0);

        let y_min = data_tol
            .iter()
            .chain(&data_seq)
            .chain(&data_cot)
            .map(|d| d.1)
            .fold(f64::INFINITY, f64::min);
        let y_max = data_tol
            .iter()
            .chain(&data_seq)
            .chain(&data_cot)
            .map(|d| d.1)
            .fold(f64::NEG_INFINITY, f64::max);
        let pad = (y_max - y_min).abs().max(1e-30) * 0.08;
        let y_lo = y_min - pad;
        let y_hi = y_max + pad;

        let datasets = vec![
            Dataset::default()
                .name("I_tol")
                .marker(symbols::Marker::Braille)
                .graph_type(GraphType::Line)
                .style(Style::default().fg(Color::Green))
                .data(&data_tol),
            Dataset::default()
                .name("I_seq")
                .marker(symbols::Marker::Braille)
                .graph_type(GraphType::Line)
                .style(Style::default().fg(Color::Yellow))
                .data(&data_seq),
            Dataset::default()
                .name("I_cot")
                .marker(symbols::Marker::Braille)
                .graph_type(GraphType::Line)
                .style(Style::default().fg(Color::Cyan))
                .data(&data_cot),
        ];

        let x_axis = Axis::default()
            .title("Vsd (V)")
            .style(Style::default().fg(Color::Gray))
            .bounds([x_min, x_max])
            .labels(vec![
                format!("{:.2}", x_min),
                format!("{:.2}", (x_min + x_max) / 2.0),
                format!("{:.2}", x_max),
            ]);

        let y_axis = Axis::default()
            .title("I (A)")
            .style(Style::default().fg(Color::Gray))
            .bounds([y_lo, y_hi])
            .labels(vec![
                format!("{:.1e}", y_lo),
                format!("{:.1e}", (y_lo + y_hi) / 2.0),
                format!("{:.1e}", y_hi),
            ]);

        let chart = Chart::new(datasets)
            .block(block)
            .x_axis(x_axis)
            .y_axis(y_axis)
            .legend_position(Some(LegendPosition::TopLeft));

        frame.render_widget(chart, area);
    } else {
        let msg = if app.compute_running {
            "Computing..."
        } else {
            "Press Enter to compute I-V curve"
        };
        frame.render_widget(
            Paragraph::new(msg)
                .block(block)
                .style(Style::default().fg(Color::DarkGray)),
            area,
        );
    }
}

fn render_stability(frame: &mut Frame, area: Rect, app: &App) {
    let title = format!(
        " Stability  N={}  l={:.1}  T={:.1}K ",
        app.n, app.lambda, app.t_kelvin
    );

    let prog_h = if app.compute_running || app.stability_progress.1 > 0 {
        3u16
    } else {
        0
    };
    let [map_area, progress_area] =
        Layout::vertical([Constraint::Fill(1), Constraint::Length(prog_h)]).areas(area);

    let block = Block::bordered().title(Line::from(title).centered());
    let has_data = app.stability_grid.iter().any(|r| !r.is_empty());

    if has_data {
        let inner = block.inner(map_area);
        frame.render_widget(block, map_area);

        let mut log_min = f64::INFINITY;
        let mut log_max = f64::NEG_INFINITY;
        for row in &app.stability_grid {
            for &val in row {
                let lv = val.abs().max(1e-30).log10();
                if lv.is_finite() {
                    log_min = log_min.min(lv);
                    log_max = log_max.max(lv);
                }
            }
        }
        if !log_min.is_finite() {
            log_min = -30.0;
        }
        if !log_max.is_finite() {
            log_max = -20.0;
        }

        let [content_area, xlabel_row] =
            Layout::vertical([Constraint::Fill(1), Constraint::Length(1)]).areas(inner);

        let [ylabel_area, heatmap_area, cbar_strip] = Layout::horizontal([
            Constraint::Length(6),
            Constraint::Fill(1),
            Constraint::Length(9),
        ])
        .areas(content_area);

        let [colorbar_area, cblabel_area] =
            Layout::horizontal([Constraint::Length(2), Constraint::Fill(1)]).areas(cbar_strip);

        frame.render_widget(
            Heatmap {
                data: &app.stability_grid,
                log_min,
                log_max,
            },
            heatmap_area,
        );

        frame.render_widget(Colorbar, colorbar_area);

        let vsd_max_val = app
            .stability_vsd_vals
            .last()
            .copied()
            .unwrap_or(app.stab_vsd_max);
        let label_style = Style::new().fg(Color::Gray);

        if ylabel_area.height >= 3 {
            let lw = ylabel_area.width.saturating_sub(1);
            let top_r = Rect::new(ylabel_area.x, ylabel_area.y, lw, 1);
            let mid_r = Rect::new(ylabel_area.x, ylabel_area.y + ylabel_area.height / 2, lw, 1);
            let bot_r = Rect::new(ylabel_area.x, ylabel_area.bottom() - 1, lw, 1);
            frame.render_widget(
                Paragraph::new(format!("{:.2}", vsd_max_val))
                    .style(label_style)
                    .alignment(Alignment::Right),
                top_r,
            );
            frame.render_widget(
                Paragraph::new(format!("{:.2}", vsd_max_val / 2.0))
                    .style(label_style)
                    .alignment(Alignment::Right),
                mid_r,
            );
            frame.render_widget(
                Paragraph::new("0.00".to_string())
                    .style(label_style)
                    .alignment(Alignment::Right),
                bot_r,
            );
        }

        let vg_lo = app.stability_vg_vals.first().copied().unwrap_or(app.vg_min);
        let vg_hi = app.stability_vg_vals.last().copied().unwrap_or(app.vg_max);
        let vg_mid = (vg_lo + vg_hi) / 2.0;
        let xl = Rect::new(heatmap_area.x, xlabel_row.y, heatmap_area.width, 1);

        if xl.width >= 20 {
            let s_left = format!("{:.2}", vg_lo);
            let s_mid = format!("{:.2}", vg_mid);
            let s_right = format!("{:.2}", vg_hi);

            let left_r = Rect::new(xl.x, xl.y, s_left.len() as u16, 1);
            frame.render_widget(Paragraph::new(s_left).style(label_style), left_r);

            let mid_w = s_mid.len() as u16;
            let mid_x = xl.x + xl.width / 2 - mid_w / 2;
            frame.render_widget(
                Paragraph::new(s_mid).style(label_style),
                Rect::new(mid_x, xl.y, mid_w, 1),
            );

            let right_w = s_right.len() as u16;
            frame.render_widget(
                Paragraph::new(s_right)
                    .style(label_style)
                    .alignment(Alignment::Right),
                Rect::new(xl.right().saturating_sub(right_w), xl.y, right_w, 1),
            );
        }

        if cblabel_area.height >= 2 {
            let cb_top = Rect::new(cblabel_area.x, cblabel_area.y, cblabel_area.width, 1);
            let cb_bot = Rect::new(
                cblabel_area.x,
                cblabel_area.bottom().saturating_sub(1),
                cblabel_area.width,
                1,
            );
            frame.render_widget(
                Paragraph::new(format!("{:.0}", log_max)).style(label_style),
                cb_top,
            );
            frame.render_widget(
                Paragraph::new(format!("{:.0}", log_min)).style(label_style),
                cb_bot,
            );
            if cblabel_area.height >= 4 {
                let cb_title = Rect::new(cblabel_area.x, cblabel_area.y + 1, cblabel_area.width, 1);
                frame.render_widget(
                    Paragraph::new("lg|I|").style(Style::new().fg(Color::DarkGray)),
                    cb_title,
                );
            }
        }
    } else {
        let msg = if app.compute_running {
            "Computing stability diagram..."
        } else {
            "Press Enter to compute stability diagram"
        };
        frame.render_widget(
            Paragraph::new(msg)
                .block(block)
                .style(Style::default().fg(Color::DarkGray)),
            map_area,
        );
    }

    if prog_h > 0 {
        let (done, total) = app.stability_progress;
        let ratio = if total > 0 {
            done as f64 / total as f64
        } else {
            0.0
        };
        let label = format!("{}/{} ({:.0}%)", done, total, ratio * 100.0);
        let gauge = Gauge::default()
            .block(Block::bordered().title("Progress"))
            .gauge_style(Style::default().fg(Color::Cyan))
            .ratio(ratio.clamp(0.0, 1.0))
            .label(label);
        frame.render_widget(gauge, progress_area);
    }
}

fn render_help(frame: &mut Frame, area: Rect, app: &App) {
    let help = if app.compute_running {
        " Esc Cancel | q Quit "
    } else {
        match app.mode {
            AppMode::IVCurve => {
                " Up/Dn Select | Lt/Rt Adjust (Shift=fine) auto-recompute | Tab Mode | e Export | q Quit "
            }
            AppMode::Stability => {
                " Up/Dn Select | Lt/Rt Adjust | Enter Run | Tab Mode | e Export | q Quit "
            }
        }
    };
    frame.render_widget(
        Paragraph::new(help).style(Style::default().fg(Color::DarkGray)),
        area,
    );
}

// ---------------------------------------------------------------------------
// Heatmap widget — half-block characters for 2x vertical resolution
// ---------------------------------------------------------------------------

struct Heatmap<'a> {
    data: &'a [Vec<f64>],
    log_min: f64,
    log_max: f64,
}

impl Widget for Heatmap<'_> {
    fn render(self, area: Rect, buf: &mut Buffer) {
        let n_cols = self.data.len();
        if n_cols == 0 {
            return;
        }
        let n_data_rows = self.data.iter().map(|r| r.len()).max().unwrap_or(0);
        if n_data_rows == 0 {
            return;
        }

        let dw = area.width as usize;
        let dh = area.height as usize;
        let pixel_h = dh * 2;
        let log_range = (self.log_max - self.log_min).max(1e-10);

        for dy in 0..dh {
            for dx in 0..dw {
                let col = if dw > 1 {
                    (dx * (n_cols - 1) / (dw - 1)).min(n_cols - 1)
                } else {
                    0
                };

                let top_pixel = dy * 2;
                let top_row = if pixel_h > 1 {
                    ((pixel_h - 1 - top_pixel) * (n_data_rows - 1) / (pixel_h - 1))
                        .min(n_data_rows - 1)
                } else {
                    0
                };

                let bot_pixel = dy * 2 + 1;
                let bot_row = if bot_pixel < pixel_h && pixel_h > 1 {
                    ((pixel_h - 1 - bot_pixel) * (n_data_rows - 1) / (pixel_h - 1))
                        .min(n_data_rows - 1)
                } else {
                    top_row
                };

                let row_data = &self.data[col];
                let top_val = if top_row < row_data.len() {
                    row_data[top_row].abs().max(1e-30).log10()
                } else {
                    self.log_min
                };
                let bot_val = if bot_row < row_data.len() {
                    row_data[bot_row].abs().max(1e-30).log10()
                } else {
                    self.log_min
                };

                let top_t = ((top_val - self.log_min) / log_range).clamp(0.0, 1.0);
                let bot_t = ((bot_val - self.log_min) / log_range).clamp(0.0, 1.0);

                let x = area.x + dx as u16;
                let y = area.y + dy as u16;
                if x < area.right() && y < area.bottom() {
                    if let Some(cell) = buf.cell_mut(ratatui::layout::Position::new(x, y)) {
                        cell.set_symbol("\u{2580}")
                            .set_style(Style::new().fg(viridis(top_t)).bg(viridis(bot_t)));
                    }
                }
            }
        }
    }
}

struct Colorbar;

impl Widget for Colorbar {
    fn render(self, area: Rect, buf: &mut Buffer) {
        let dh = area.height as usize;
        if dh == 0 {
            return;
        }
        let pixel_h = dh * 2;

        for dy in 0..dh {
            let top_pixel = dy * 2;
            let top_t = if pixel_h > 1 {
                (pixel_h - 1 - top_pixel) as f64 / (pixel_h - 1) as f64
            } else {
                0.5
            };

            let bot_pixel = dy * 2 + 1;
            let bot_t = if bot_pixel < pixel_h && pixel_h > 1 {
                (pixel_h - 1 - bot_pixel) as f64 / (pixel_h - 1) as f64
            } else {
                top_t
            };

            for dx in 0..area.width {
                let x = area.x + dx;
                let y = area.y + dy as u16;
                if let Some(cell) = buf.cell_mut(ratatui::layout::Position::new(x, y)) {
                    cell.set_symbol("\u{2580}")
                        .set_style(Style::new().fg(viridis(top_t)).bg(viridis(bot_t)));
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Event handling
// ---------------------------------------------------------------------------

fn handle_key(app: &mut App, key: event::KeyEvent) {
    match key.code {
        KeyCode::Char('q') | KeyCode::Char('Q') => {
            app.should_quit = true;
        }
        KeyCode::Esc => {
            if app.compute_running {
                app.cancel_compute();
            } else {
                app.should_quit = true;
            }
        }
        KeyCode::Tab => {
            if !app.compute_running {
                app.mode = match app.mode {
                    AppMode::IVCurve => AppMode::Stability,
                    AppMode::Stability => AppMode::IVCurve,
                };
                app.selected_param = 0;
                app.iv_recompute_at = None;
            }
        }
        KeyCode::Up => {
            if app.selected_param > 0 {
                app.selected_param -= 1;
            }
        }
        KeyCode::Down => {
            let max = app.visible_params().len().saturating_sub(1);
            if app.selected_param < max {
                app.selected_param += 1;
            }
        }
        KeyCode::Right | KeyCode::Left => {
            if app.mode == AppMode::Stability && app.compute_running {
                return;
            }
            let params = app.visible_params();
            if let Some(&id) = params.get(app.selected_param) {
                let dir = if key.code == KeyCode::Right { 1 } else { -1 };
                let fine = key.modifiers.contains(KeyModifiers::SHIFT);
                app.adjust_param(id, dir, fine);
                app.status_msg = None;
                if app.mode == AppMode::IVCurve {
                    app.iv_recompute_at = Some(Instant::now() + Duration::from_millis(300));
                }
            }
        }
        KeyCode::Enter => {
            app.iv_recompute_at = None;
            if app.compute_running {
                app.cancel_compute();
            } else {
                match app.mode {
                    AppMode::IVCurve => app.start_iv(),
                    AppMode::Stability => app.start_stability(),
                }
            }
        }
        KeyCode::Char('e') | KeyCode::Char('E') => {
            if !app.compute_running {
                let result = match app.mode {
                    AppMode::Stability => app.export_stability(),
                    AppMode::IVCurve => app.export_iv(),
                };
                match result {
                    Ok(path) => app.status_msg = Some((format!("Exported: {}", path), false)),
                    Err(e) => app.status_msg = Some((format!("Export failed: {}", e), true)),
                }
            }
        }
        _ => {}
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
        terminal.draw(|frame| render(frame, &app))?;

        let timeout = TICK_RATE.saturating_sub(last_tick.elapsed());
        if event::poll(timeout)? {
            if let Event::Key(key) = event::read()? {
                if key.kind == KeyEventKind::Press {
                    handle_key(&mut app, key);
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
