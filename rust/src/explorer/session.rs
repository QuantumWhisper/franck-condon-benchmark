use std::time::Duration;

use serde::{Deserialize, Serialize};

use franck_condon::simulate::SimulationResult;

use super::{App, AppMode, DisplayMode};

pub(crate) const DEFAULT_SESSION_PATH: &str = "fc_session.json";

const SESSION_VERSION: u32 = 1;

// ---------------------------------------------------------------------------
// Custom serde for f64 values that may be Infinity (JSON doesn't support it)
// Matches the project convention: tau="Inf" in benchmark spec JSON files.
// ---------------------------------------------------------------------------

mod f64_inf {
    use serde::{self, Deserialize, Deserializer, Serializer};

    pub fn serialize<S>(val: &f64, ser: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        if val.is_infinite() && val.is_sign_positive() {
            ser.serialize_str("Inf")
        } else if val.is_infinite() {
            ser.serialize_str("-Inf")
        } else if val.is_nan() {
            ser.serialize_str("NaN")
        } else {
            ser.serialize_f64(*val)
        }
    }

    pub fn deserialize<'de, D>(de: D) -> Result<f64, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(Deserialize)]
        #[serde(untagged)]
        enum FloatOrStr {
            Float(f64),
            Str(String),
        }
        match FloatOrStr::deserialize(de)? {
            FloatOrStr::Float(v) => Ok(v),
            FloatOrStr::Str(s) => match s.as_str() {
                "Inf" | "inf" | "Infinity" | "infinity" => Ok(f64::INFINITY),
                "-Inf" | "-inf" | "-Infinity" | "-infinity" => Ok(f64::NEG_INFINITY),
                "NaN" | "nan" => Ok(f64::NAN),
                other => other.parse::<f64>().map_err(serde::de::Error::custom),
            },
        }
    }
}

// Duration stored as f64 seconds since std::time::Duration isn't JSON-friendly
#[derive(Serialize, Deserialize)]
struct SavedIVData {
    result: SimulationResult,
    elapsed_secs: f64,
}

#[derive(Serialize, Deserialize)]
pub(crate) struct Session {
    version: u32,

    // Current mode and display
    mode: AppMode,
    display_mode: DisplayMode,

    // Physics parameters
    n: usize,
    lambda: f64,
    t_kelvin: f64,
    vmode: f64,
    alpha_l: f64,
    alpha_r: f64,
    eta: f64,
    vg: f64,
    #[serde(with = "f64_inf")]
    tau: f64,

    // I-V sweep settings
    vsd_min: f64,
    vsd_max: f64,
    vsd_step: f64,

    // Stability grid settings
    stab_vsd_max: f64,
    n_vsd: usize,
    vg_min: f64,
    vg_max: f64,
    n_vg: usize,

    // Temperature diagram settings
    temp_vsd_max: f64,
    temp_n_vsd: usize,
    t_min: f64,
    t_max: f64,
    n_t: usize,

    // Computed I-V data
    iv_data: Option<SavedIVData>,

    // Computed stability diagram
    stability_grid: Vec<Vec<f64>>,
    stability_vg_vals: Vec<f64>,
    stability_vsd_vals: Vec<f64>,
    stability_elapsed_secs: Option<f64>,

    // Computed temperature diagram
    temperature_grid: Vec<Vec<f64>>,
    temperature_t_vals: Vec<f64>,
    temperature_vsd_vals: Vec<f64>,
    temperature_elapsed_secs: Option<f64>,

    // Display settings
    contrast_min: Option<f64>,
    contrast_max: Option<f64>,
    gamma: f64,
}

impl Session {
    pub(crate) fn from_app(app: &App) -> Self {
        Session {
            version: SESSION_VERSION,
            mode: app.mode,
            display_mode: app.display_mode,
            n: app.n,
            lambda: app.lambda,
            t_kelvin: app.t_kelvin,
            vmode: app.vmode,
            alpha_l: app.alpha_l,
            alpha_r: app.alpha_r,
            eta: app.eta,
            vg: app.vg,
            tau: app.tau,
            vsd_min: app.vsd_min,
            vsd_max: app.vsd_max,
            vsd_step: app.vsd_step,
            stab_vsd_max: app.stab_vsd_max,
            n_vsd: app.n_vsd,
            vg_min: app.vg_min,
            vg_max: app.vg_max,
            n_vg: app.n_vg,
            temp_vsd_max: app.temp_vsd_max,
            temp_n_vsd: app.temp_n_vsd,
            t_min: app.t_min,
            t_max: app.t_max,
            n_t: app.n_t,
            iv_data: app.iv_data.as_ref().map(|(result, elapsed)| SavedIVData {
                result: result.clone(),
                elapsed_secs: elapsed.as_secs_f64(),
            }),
            stability_grid: app.stability_grid.clone(),
            stability_vg_vals: app.stability_vg_vals.clone(),
            stability_vsd_vals: app.stability_vsd_vals.clone(),
            stability_elapsed_secs: app.stability_elapsed.map(|d| d.as_secs_f64()),
            temperature_grid: app.temperature_grid.clone(),
            temperature_t_vals: app.temperature_t_vals.clone(),
            temperature_vsd_vals: app.temperature_vsd_vals.clone(),
            temperature_elapsed_secs: app.temperature_elapsed.map(|d| d.as_secs_f64()),
            contrast_min: app.contrast_min,
            contrast_max: app.contrast_max,
            gamma: app.gamma,
        }
    }

    pub(crate) fn apply_to(self, app: &mut App) {
        app.mode = self.mode;
        app.display_mode = self.display_mode;
        app.n = self.n;
        app.lambda = self.lambda;
        app.t_kelvin = self.t_kelvin;
        app.vmode = self.vmode;
        app.alpha_l = self.alpha_l;
        app.alpha_r = self.alpha_r;
        app.eta = self.eta;
        app.vg = self.vg;
        app.tau = self.tau;
        app.vsd_min = self.vsd_min;
        app.vsd_max = self.vsd_max;
        app.vsd_step = self.vsd_step;
        app.stab_vsd_max = self.stab_vsd_max;
        app.n_vsd = self.n_vsd;
        app.vg_min = self.vg_min;
        app.vg_max = self.vg_max;
        app.n_vg = self.n_vg;
        app.temp_vsd_max = self.temp_vsd_max;
        app.temp_n_vsd = self.temp_n_vsd;
        app.t_min = self.t_min;
        app.t_max = self.t_max;
        app.n_t = self.n_t;

        app.iv_data = self
            .iv_data
            .map(|saved| (saved.result, Duration::from_secs_f64(saved.elapsed_secs)));

        app.stability_grid = self.stability_grid;
        app.stability_vg_vals = self.stability_vg_vals;
        app.stability_vsd_vals = self.stability_vsd_vals;
        app.stability_elapsed = self.stability_elapsed_secs.map(Duration::from_secs_f64);
        let stab_done = app.stability_grid.iter().filter(|r| !r.is_empty()).count();
        app.stability_progress = if app.stability_grid.is_empty() {
            (0, 0)
        } else {
            (stab_done, app.stability_grid.len())
        };

        app.temperature_grid = self.temperature_grid;
        app.temperature_t_vals = self.temperature_t_vals;
        app.temperature_vsd_vals = self.temperature_vsd_vals;
        app.temperature_elapsed = self.temperature_elapsed_secs.map(Duration::from_secs_f64);
        let temp_done = app
            .temperature_grid
            .iter()
            .filter(|r| !r.is_empty())
            .count();
        app.temperature_progress = if app.temperature_grid.is_empty() {
            (0, 0)
        } else {
            (temp_done, app.temperature_grid.len())
        };

        app.contrast_min = self.contrast_min;
        app.contrast_max = self.contrast_max;
        app.gamma = self.gamma;

        // Reset ephemeral UI state
        app.selected_param = 0;
        app.cursor_mode = false;
        app.show_line_profile = false;
        app.iv_recompute_at = None;
    }
}

// ---------------------------------------------------------------------------
// Public save / load API
// ---------------------------------------------------------------------------

pub(crate) fn save_session(app: &App, path: &str) -> std::io::Result<String> {
    let session = Session::from_app(app);
    let json = serde_json::to_string_pretty(&session)
        .map_err(|e| std::io::Error::other(format!("Serialize: {}", e)))?;
    std::fs::write(path, json)?;
    Ok(path.to_string())
}

pub(crate) fn load_session(path: &str) -> std::io::Result<Session> {
    let json = std::fs::read_to_string(path)?;
    let session: Session = serde_json::from_str(&json)
        .map_err(|e| std::io::Error::other(format!("Deserialize: {}", e)))?;
    if session.version > SESSION_VERSION {
        return Err(std::io::Error::other(format!(
            "Session version {} is newer than supported version {}",
            session.version, SESSION_VERSION
        )));
    }
    Ok(session)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_f64_inf_roundtrip() {
        #[derive(Serialize, Deserialize, Debug)]
        struct T {
            #[serde(with = "f64_inf")]
            val: f64,
        }
        for &v in &[1.0, -3.14, 0.0, f64::INFINITY, f64::NEG_INFINITY] {
            let json = serde_json::to_string(&T { val: v }).unwrap();
            let back: T = serde_json::from_str(&json).unwrap();
            if v.is_infinite() {
                assert_eq!(v.is_sign_positive(), back.val.is_sign_positive());
                assert!(back.val.is_infinite());
            } else {
                assert!((v - back.val).abs() < 1e-15);
            }
        }
    }

    #[test]
    fn test_f64_inf_deserialize_strings() {
        #[derive(Deserialize)]
        struct T {
            #[serde(with = "f64_inf")]
            val: f64,
        }
        for s in &["\"Inf\"", "\"inf\"", "\"Infinity\""] {
            let json = format!("{{\"val\": {}}}", s);
            let t: T = serde_json::from_str(&json).unwrap();
            assert!(t.val.is_infinite() && t.val.is_sign_positive());
        }
        for s in &["\"-Inf\"", "\"-inf\""] {
            let json = format!("{{\"val\": {}}}", s);
            let t: T = serde_json::from_str(&json).unwrap();
            assert!(t.val.is_infinite() && t.val.is_sign_negative());
        }
    }

    #[test]
    fn test_session_version_reject_newer() {
        let json = r#"{"version": 999, "mode": "IVCurve", "display_mode": "Current",
            "n": 6, "lambda": 5.0, "t_kelvin": 4.2, "vmode": 0.073,
            "alpha_l": 0.02, "alpha_r": 0.02, "eta": 0.5, "vg": 0.0, "tau": "Inf",
            "vsd_min": 0.0, "vsd_max": 0.6, "vsd_step": 0.003,
            "stab_vsd_max": 0.6, "n_vsd": 101, "vg_min": -0.3, "vg_max": 0.3, "n_vg": 61,
            "temp_vsd_max": 0.6, "temp_n_vsd": 101, "t_min": 1.0, "t_max": 50.0, "n_t": 50,
            "iv_data": null,
            "stability_grid": [], "stability_vg_vals": [], "stability_vsd_vals": [],
            "stability_elapsed_secs": null,
            "temperature_grid": [], "temperature_t_vals": [], "temperature_vsd_vals": [],
            "temperature_elapsed_secs": null,
            "contrast_min": null, "contrast_max": null, "gamma": 1.0}"#;

        let tmp = std::env::temp_dir().join("fc_test_version.json");
        std::fs::write(&tmp, json).unwrap();
        let result = load_session(tmp.to_str().unwrap());
        assert!(result.is_err());
        let err = result.err().unwrap();
        assert!(err.to_string().contains("newer"), "got: {}", err);
        let _ = std::fs::remove_file(&tmp);
    }

    #[test]
    fn test_session_roundtrip_minimal() {
        let session = Session {
            version: SESSION_VERSION,
            mode: AppMode::IVCurve,
            display_mode: DisplayMode::Current,
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
            temp_vsd_max: 0.6,
            temp_n_vsd: 101,
            t_min: 1.0,
            t_max: 50.0,
            n_t: 50,
            iv_data: None,
            stability_grid: vec![],
            stability_vg_vals: vec![],
            stability_vsd_vals: vec![],
            stability_elapsed_secs: None,
            temperature_grid: vec![],
            temperature_t_vals: vec![],
            temperature_vsd_vals: vec![],
            temperature_elapsed_secs: None,
            contrast_min: None,
            contrast_max: None,
            gamma: 1.0,
        };

        let json = serde_json::to_string_pretty(&session).unwrap();
        let back: Session = serde_json::from_str(&json).unwrap();

        assert_eq!(back.n, 6);
        assert!(back.tau.is_infinite());
        assert!(back.iv_data.is_none());
        assert_eq!(back.mode, AppMode::IVCurve);
        assert_eq!(back.display_mode, DisplayMode::Current);
    }

    #[test]
    fn test_session_roundtrip_with_data() {
        let session = Session {
            version: SESSION_VERSION,
            mode: AppMode::Stability,
            display_mode: DisplayMode::Conductance,
            n: 6,
            lambda: 5.0,
            t_kelvin: 4.2,
            vmode: 0.073,
            alpha_l: 0.02,
            alpha_r: 0.02,
            eta: 0.5,
            vg: 0.0,
            tau: 1e-9,
            vsd_min: 0.0,
            vsd_max: 0.6,
            vsd_step: 0.003,
            stab_vsd_max: 0.6,
            n_vsd: 3,
            vg_min: -0.3,
            vg_max: 0.3,
            n_vg: 2,
            temp_vsd_max: 0.6,
            temp_n_vsd: 3,
            t_min: 1.0,
            t_max: 50.0,
            n_t: 2,
            iv_data: Some(SavedIVData {
                result: SimulationResult {
                    vsd: vec![0.0, 0.3, 0.6],
                    i_tol: vec![1e-10, 2e-8, 5e-7],
                    i_seq: vec![0.5e-10, 1e-8, 3e-7],
                    i_cot: vec![0.5e-10, 1e-8, 2e-7],
                },
                elapsed_secs: 1.234,
            }),
            stability_grid: vec![vec![1e-10, 2e-9, 3e-8], vec![4e-10, 5e-9, 6e-8]],
            stability_vg_vals: vec![-0.3, 0.3],
            stability_vsd_vals: vec![0.0, 0.3, 0.6],
            stability_elapsed_secs: Some(42.5),
            temperature_grid: vec![vec![1e-11, 2e-10], vec![3e-11, 4e-10]],
            temperature_t_vals: vec![1.0, 50.0],
            temperature_vsd_vals: vec![0.0, 0.6],
            temperature_elapsed_secs: Some(15.3),
            contrast_min: Some(-25.0),
            contrast_max: Some(-5.0),
            gamma: 0.8,
        };

        let json = serde_json::to_string_pretty(&session).unwrap();
        let back: Session = serde_json::from_str(&json).unwrap();

        assert_eq!(back.mode, AppMode::Stability);
        assert_eq!(back.display_mode, DisplayMode::Conductance);
        assert!((back.tau - 1e-9).abs() < 1e-20);
        assert!(back.iv_data.is_some());
        let iv = back.iv_data.unwrap();
        assert_eq!(iv.result.vsd.len(), 3);
        assert!((iv.elapsed_secs - 1.234).abs() < 1e-10);
        assert_eq!(back.stability_grid.len(), 2);
        assert_eq!(back.stability_grid[0].len(), 3);
        assert!((back.gamma - 0.8).abs() < 1e-10);
        assert_eq!(back.contrast_min, Some(-25.0));
    }

    #[test]
    fn test_save_and_load_file() {
        let tmp = std::env::temp_dir().join("fc_test_save_load.json");
        let path = tmp.to_str().unwrap();

        let session = Session {
            version: SESSION_VERSION,
            mode: AppMode::Temperature,
            display_mode: DisplayMode::Iets,
            n: 10,
            lambda: 3.0,
            t_kelvin: 10.0,
            vmode: 0.073,
            alpha_l: 0.02,
            alpha_r: 0.02,
            eta: 0.5,
            vg: -0.1,
            tau: f64::INFINITY,
            vsd_min: 0.0,
            vsd_max: 1.0,
            vsd_step: 0.005,
            stab_vsd_max: 1.0,
            n_vsd: 201,
            vg_min: -0.5,
            vg_max: 0.5,
            n_vg: 101,
            temp_vsd_max: 1.0,
            temp_n_vsd: 201,
            t_min: 2.0,
            t_max: 100.0,
            n_t: 99,
            iv_data: None,
            stability_grid: vec![],
            stability_vg_vals: vec![],
            stability_vsd_vals: vec![],
            stability_elapsed_secs: None,
            temperature_grid: vec![],
            temperature_t_vals: vec![],
            temperature_vsd_vals: vec![],
            temperature_elapsed_secs: None,
            contrast_min: None,
            contrast_max: None,
            gamma: 1.0,
        };

        let json = serde_json::to_string_pretty(&session).unwrap();
        std::fs::write(path, &json).unwrap();

        let back = load_session(path).unwrap();
        assert_eq!(back.mode, AppMode::Temperature);
        assert_eq!(back.n, 10);
        assert!(back.tau.is_infinite());

        let _ = std::fs::remove_file(&tmp);
    }
}
