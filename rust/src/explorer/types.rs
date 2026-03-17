use std::time::Duration;

use franck_condon::simulate::SimulationResult;

use super::derivatives::*;

#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum AppMode {
    IVCurve,
    Stability,
    Temperature,
}

impl AppMode {
    pub(crate) fn label(self) -> &'static str {
        match self {
            Self::IVCurve => "I-V Curve",
            Self::Stability => "Stability Diagram",
            Self::Temperature => "Temperature Diagram",
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum DisplayMode {
    Current,
    Conductance,
    Iets,
    NormalizedIets,
}

impl DisplayMode {
    pub(crate) fn next(self) -> Self {
        match self {
            Self::Current => Self::Conductance,
            Self::Conductance => Self::Iets,
            Self::Iets => Self::NormalizedIets,
            Self::NormalizedIets => Self::Current,
        }
    }

    pub(crate) fn suffix(self) -> &'static str {
        match self {
            Self::Current => "",
            Self::Conductance => " [dI/dV]",
            Self::Iets => " [IETS]",
            Self::NormalizedIets => " [nIETS]",
        }
    }

    pub(crate) fn y_label(self) -> &'static str {
        match self {
            Self::Current => "I (A)",
            Self::Conductance => "dI/dV (S)",
            Self::Iets => "d\u{00b2}I/dV\u{00b2} (S/V)",
            Self::NormalizedIets => "IETS/G (1/V)",
        }
    }

    pub(crate) fn dataset_names(self) -> (&'static str, &'static str, &'static str) {
        match self {
            Self::Current => ("I_tol", "I_seq", "I_cot"),
            Self::Conductance => ("G_tol", "G_seq", "G_cot"),
            Self::Iets => ("d\u{00b2}_tol", "d\u{00b2}_seq", "d\u{00b2}_cot"),
            Self::NormalizedIets => ("nI_tol", "nI_seq", "nI_cot"),
        }
    }

    pub(crate) fn colorbar_label(self) -> &'static str {
        match self {
            Self::Current => "lg|I|",
            Self::Conductance => "lg|G|",
            Self::Iets => "lg|d\u{00b2}I|",
            Self::NormalizedIets => "lg|nI|",
        }
    }

    pub(crate) fn iv_csv_header(self) -> &'static str {
        match self {
            Self::Current => "Vsd_V,I_tol_A,I_seq_A,I_cot_A",
            Self::Conductance => "Vsd_V,dIdV_tol_S,dIdV_seq_S,dIdV_cot_S",
            Self::Iets => "Vsd_V,d2IdV2_tol,d2IdV2_seq,d2IdV2_cot",
            Self::NormalizedIets => "Vsd_V,nIETS_tol,nIETS_seq,nIETS_cot",
        }
    }

    pub(crate) fn grid_csv_column(self) -> &'static str {
        match self {
            Self::Current => "I_tol_A",
            Self::Conductance => "dIdV_S",
            Self::Iets => "d2IdV2",
            Self::NormalizedIets => "nIETS",
        }
    }

    pub(crate) fn iv_export_filename(self) -> &'static str {
        match self {
            Self::Current => "iv_export.csv",
            Self::Conductance => "iv_didv_export.csv",
            Self::Iets => "iv_iets_export.csv",
            Self::NormalizedIets => "iv_niets_export.csv",
        }
    }

    pub(crate) fn stability_export_filename(self) -> &'static str {
        match self {
            Self::Current => "stability_export.csv",
            Self::Conductance => "stability_didv_export.csv",
            Self::Iets => "stability_iets_export.csv",
            Self::NormalizedIets => "stability_niets_export.csv",
        }
    }

    pub(crate) fn temperature_export_filename(self) -> &'static str {
        match self {
            Self::Current => "temperature_export.csv",
            Self::Conductance => "temperature_didv_export.csv",
            Self::Iets => "temperature_iets_export.csv",
            Self::NormalizedIets => "temperature_niets_export.csv",
        }
    }

    pub(crate) fn transform(self, vsd: &[f64], current: &[f64]) -> Vec<f64> {
        match self {
            Self::Current => current.to_vec(),
            Self::Conductance => compute_didv(vsd, current),
            Self::Iets => compute_d2idv2(vsd, current),
            Self::NormalizedIets => compute_normalized_iets(vsd, current),
        }
    }

    pub(crate) fn transform_grid(self, grid: &[Vec<f64>], vsd_vals: &[f64]) -> Vec<Vec<f64>> {
        match self {
            Self::Current => grid.to_vec(),
            Self::Conductance => compute_didv_grid(grid, vsd_vals),
            Self::Iets => compute_d2idv2_grid(grid, vsd_vals),
            Self::NormalizedIets => compute_normalized_iets_grid(grid, vsd_vals),
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum ParamId {
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
    TempVsdMax,
    TempNVsd,
    TMin,
    TMax,
    NT,
}

pub(crate) enum ComputeMsg {
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
    TemperatureRow {
        t_idx: usize,
        i_tol: Vec<f64>,
    },
    TemperatureDone {
        elapsed: Duration,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_display_mode_cycle() {
        let mut dm = DisplayMode::Current;
        dm = dm.next();
        assert_eq!(dm, DisplayMode::Conductance);
        dm = dm.next();
        assert_eq!(dm, DisplayMode::Iets);
        dm = dm.next();
        assert_eq!(dm, DisplayMode::NormalizedIets);
        dm = dm.next();
        assert_eq!(dm, DisplayMode::Current);
    }

    #[test]
    fn test_app_mode_labels() {
        assert!(!AppMode::IVCurve.label().is_empty());
        assert!(!AppMode::Stability.label().is_empty());
        assert!(!AppMode::Temperature.label().is_empty());
    }

    #[test]
    fn test_display_mode_suffix() {
        assert_eq!(DisplayMode::Current.suffix(), "");
        assert!(DisplayMode::Conductance.suffix().contains("dI/dV"));
    }

    #[test]
    fn test_colorbar_labels() {
        assert!(DisplayMode::Current.colorbar_label().contains("I"));
        assert!(DisplayMode::Conductance.colorbar_label().contains("G"));
    }

    #[test]
    fn test_export_filenames_non_empty() {
        assert!(!DisplayMode::Iets.iv_export_filename().is_empty());
        assert!(!DisplayMode::Iets.stability_export_filename().is_empty());
        assert!(!DisplayMode::Iets.temperature_export_filename().is_empty());
    }
}
