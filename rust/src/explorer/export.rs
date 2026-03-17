use std::io::{self, Write};
use std::process::{Command, Stdio};

use super::widgets::VIRIDIS;
use super::{App, AppMode, DisplayMode};

pub(crate) fn export_stability(app: &App) -> io::Result<String> {
    let buf: Vec<Vec<f64>>;
    let grid: &[Vec<f64>] = if app.display_mode == super::DisplayMode::Current {
        &app.stability_grid
    } else {
        buf = app
            .display_mode
            .transform_grid(&app.stability_grid, &app.stability_vsd_vals);
        &buf
    };
    let filename = app.display_mode.stability_export_filename();
    let header = format!("Vg_V,Vsd_V,{}", app.display_mode.grid_csv_column());
    let mut f = std::fs::File::create(filename)?;
    writeln!(f, "{}", header)?;
    for (i, &vg) in app.stability_vg_vals.iter().enumerate() {
        if i >= grid.len() || grid[i].is_empty() {
            continue;
        }
        for (j, &vsd) in app.stability_vsd_vals.iter().enumerate() {
            if j < grid[i].len() {
                writeln!(f, "{:.6e},{:.6e},{:.6e}", vg, vsd, grid[i][j])?;
            }
        }
    }
    Ok(filename.into())
}

pub(crate) fn export_iv(app: &App) -> io::Result<String> {
    if let Some((ref data, _)) = app.iv_data {
        let filename = app.display_mode.iv_export_filename();
        let v1 = app.display_mode.transform(&data.vsd, &data.i_tol);
        let v2 = app.display_mode.transform(&data.vsd, &data.i_seq);
        let v3 = app.display_mode.transform(&data.vsd, &data.i_cot);
        let mut f = std::fs::File::create(filename)?;
        writeln!(f, "{}", app.display_mode.iv_csv_header())?;
        for i in 0..data.vsd.len() {
            writeln!(
                f,
                "{:.6e},{:.6e},{:.6e},{:.6e}",
                data.vsd[i], v1[i], v2[i], v3[i]
            )?;
        }
        Ok(filename.into())
    } else {
        Err(io::Error::other("No I-V data"))
    }
}

pub(crate) fn export_temperature(app: &App) -> io::Result<String> {
    let buf: Vec<Vec<f64>>;
    let grid: &[Vec<f64>] = if app.display_mode == super::DisplayMode::Current {
        &app.temperature_grid
    } else {
        buf = app
            .display_mode
            .transform_grid(&app.temperature_grid, &app.temperature_vsd_vals);
        &buf
    };
    let filename = app.display_mode.temperature_export_filename();
    let header = format!("T_K,Vsd_V,{}", app.display_mode.grid_csv_column());
    let mut f = std::fs::File::create(filename)?;
    writeln!(f, "{}", header)?;
    for (i, &t) in app.temperature_t_vals.iter().enumerate() {
        if i >= grid.len() || grid[i].is_empty() {
            continue;
        }
        for (j, &vsd) in app.temperature_vsd_vals.iter().enumerate() {
            if j < grid[i].len() {
                writeln!(f, "{:.6e},{:.6e},{:.6e}", t, vsd, grid[i][j])?;
            }
        }
    }
    Ok(filename.into())
}

fn display_mode_tag(display_mode: DisplayMode) -> &'static str {
    match display_mode {
        DisplayMode::Current => "current",
        DisplayMode::Conductance => "conductance",
        DisplayMode::Iets => "iets",
        DisplayMode::NormalizedIets => "niets",
    }
}

fn iv_plot_paths(display_mode: DisplayMode) -> (String, String) {
    let tag = display_mode_tag(display_mode);
    (format!("iv_{}.pdf", tag), format!("iv_{}.png", tag))
}

fn heatmap_plot_paths(mode: AppMode, display_mode: DisplayMode) -> (String, String) {
    let mode_tag = match mode {
        AppMode::Stability => "stability",
        AppMode::Temperature => "temperature",
        AppMode::IVCurve => "iv",
    };
    let dm_tag = display_mode_tag(display_mode);
    (
        format!("{}_{}.pdf", mode_tag, dm_tag),
        format!("{}_{}.png", mode_tag, dm_tag),
    )
}

fn viridis_palette_string() -> String {
    let mut parts = Vec::new();
    let n = VIRIDIS.len();
    for (i, &(r, g, b)) in VIRIDIS.iter().enumerate() {
        let t = i as f64 / (n - 1) as f64;
        let ri = (r * 255.0) as u8;
        let gi = (g * 255.0) as u8;
        let bi = (b * 255.0) as u8;
        parts.push(format!("{:.3} '#{:02x}{:02x}{:02x}'", t, ri, gi, bi));
    }
    format!("set palette defined ({})", parts.join(", "))
}

fn run_gnuplot(script: &str) -> io::Result<()> {
    let mut child = Command::new("gnuplot")
        .stdin(Stdio::piped())
        .stdout(Stdio::null())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|_| io::Error::other("gnuplot not available"))?;

    if let Some(stdin) = child.stdin.as_mut() {
        stdin.write_all(script.as_bytes())?;
    }

    let output = child.wait_with_output()?;
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        let msg = if stderr.trim().is_empty() {
            "gnuplot failed".to_string()
        } else {
            format!("gnuplot failed: {}", stderr.trim())
        };
        return Err(io::Error::other(msg));
    }
    Ok(())
}

pub(crate) fn export_iv_gnuplot(app: &App) -> io::Result<String> {
    let (data, _) = app
        .iv_data
        .as_ref()
        .ok_or_else(|| io::Error::other("No I-V data"))?;

    let v1 = app.display_mode.transform(&data.vsd, &data.i_tol);
    let v2 = app.display_mode.transform(&data.vsd, &data.i_seq);
    let v3 = app.display_mode.transform(&data.vsd, &data.i_cot);

    let csv_path = "_gnuplot_iv_temp.csv";
    {
        let mut f = std::fs::File::create(csv_path)?;
        writeln!(f, "Vsd,total,seq,cot")?;
        for i in 0..data.vsd.len() {
            writeln!(
                f,
                "{:.6e},{:.6e},{:.6e},{:.6e}",
                data.vsd[i], v1[i], v2[i], v3[i]
            )?;
        }
    }

    let (pdf_path, png_path) = iv_plot_paths(app.display_mode);
    let (name_tol, name_seq, name_cot) = app.display_mode.dataset_names();
    let y_label = app.display_mode.y_label();
    let title = format!(
        "N={}, lambda={:.1}, T={:.1} K, hbar*omega={:.0} meV",
        app.n,
        app.lambda,
        app.t_kelvin,
        app.vmode * 1e3
    );
    let annotation = format!(
        "alpha_L={:.2}, alpha_R={:.2}, eta={:.1}, Vg={:.1} V",
        app.alpha_l, app.alpha_r, app.eta, app.vg
    );

    let script = format!(
        r#"set datafile separator ','
set key left top box opaque
set grid lw 0.5 lc rgb '#cccccc'
set xlabel 'V_{{sd}} (V)'
set ylabel '{y_label}'
set title '{title}'
set label 1 '{annotation}' at graph 0.98,0.03 right front font ',8'
set style line 1 lc rgb '#111111' pt 7 ps 0.35 lw 1.0
set style line 2 lc rgb '#0059b3' lw 1.5
set style line 3 lc rgb '#cc1f1f' lw 1.5 dt 2
set terminal pdfcairo enhanced font 'Times,10' size 12cm,9cm
set output '{pdf_path}'
plot '{csv_path}' every ::1 using 1:2 with points ls 1 title '{name_tol}', \
     '' every ::1 using 1:3 with lines ls 2 title '{name_seq}', \
     '' every ::1 using 1:4 with lines ls 3 title '{name_cot}'
set terminal pngcairo enhanced font 'Times,10' size 12cm,9cm
set output '{png_path}'
replot
unset output
"#
    );

    let run_result = run_gnuplot(&script);
    let _ = std::fs::remove_file(csv_path);
    run_result.map(|_| format!("{}, {}", pdf_path, png_path))
}

pub(crate) fn export_heatmap_gnuplot(app: &App) -> io::Result<String> {
    let (raw_grid, x_vals, y_vals, x_label) = match app.mode {
        AppMode::Stability => (
            &app.stability_grid,
            &app.stability_vg_vals,
            &app.stability_vsd_vals,
            "V_g (V)",
        ),
        AppMode::Temperature => (
            &app.temperature_grid,
            &app.temperature_t_vals,
            &app.temperature_vsd_vals,
            "T (K)",
        ),
        AppMode::IVCurve => {
            return Err(io::Error::other(
                "Heatmap export only in Stability/Temperature mode",
            ))
        }
    };

    if raw_grid.is_empty() || raw_grid.iter().all(|r| r.is_empty()) {
        return Err(io::Error::other("No heatmap data computed yet"));
    }

    let grid_buf: Vec<Vec<f64>>;
    let display_grid: &[Vec<f64>] = if app.display_mode == DisplayMode::Current {
        raw_grid
    } else {
        grid_buf = app.display_mode.transform_grid(raw_grid, y_vals);
        &grid_buf
    };

    let mut auto_min = f64::INFINITY;
    let mut auto_max = f64::NEG_INFINITY;
    for row in display_grid {
        for &val in row {
            let lv = val.abs().max(1e-30).log10();
            if lv.is_finite() {
                auto_min = auto_min.min(lv);
                auto_max = auto_max.max(lv);
            }
        }
    }
    if !auto_min.is_finite() {
        auto_min = -30.0;
    }
    if !auto_max.is_finite() {
        auto_max = -20.0;
    }

    let effective_min = app.contrast_min.unwrap_or(auto_min);
    let effective_max = app.contrast_max.unwrap_or(auto_max);

    let dat_path = "_gnuplot_heatmap_temp.dat";
    {
        let mut f = std::fs::File::create(dat_path)?;
        write!(f, "{}", display_grid.len())?;
        for &xv in x_vals {
            write!(f, " {:.6e}", xv)?;
        }
        writeln!(f)?;

        let n_rows = display_grid.iter().map(|c| c.len()).max().unwrap_or(0);
        for j in 0..n_rows {
            let yv = y_vals.get(j).copied().unwrap_or(0.0);
            write!(f, "{:.6e}", yv)?;
            for col in display_grid {
                let raw_val = col.get(j).copied().unwrap_or(0.0);
                let log_val = raw_val.abs().max(1e-30).log10();
                let range = (effective_max - effective_min).max(1e-10);
                let t = ((log_val - effective_min) / range).clamp(0.0, 1.0);
                let t_gamma = t.powf(app.gamma);
                let z = effective_min + t_gamma * range;
                write!(f, " {:.6e}", z)?;
            }
            writeln!(f)?;
        }
    }

    let (pdf_path, png_path) = heatmap_plot_paths(app.mode, app.display_mode);
    let y_label = "V_{sd} (V)";
    let cb_label = app.display_mode.colorbar_label();
    let title = match app.mode {
        AppMode::Stability => format!(
            "N={}, lambda={:.1}, T={:.1} K",
            app.n, app.lambda, app.t_kelvin
        ),
        AppMode::Temperature => {
            format!("N={}, lambda={:.1}, Vg={:.3} V", app.n, app.lambda, app.vg)
        }
        AppMode::IVCurve => String::new(),
    };

    let palette = viridis_palette_string();
    let script = format!(
        r#"set terminal pdfcairo enhanced font 'Times,10' size 12cm,9cm
set output '{pdf_path}'
{palette}
set pm3d map
set xlabel '{x_label}'
set ylabel '{y_label}'
set cblabel '{cb_label}'
set cbrange [{effective_min:.2}:{effective_max:.2}]
set title '{title}'
set grid front lw 0.3 lc rgb '#666666'
set datafile separator whitespace
splot '{dat_path}' matrix rowheaders columnheaders with image notitle
set terminal pngcairo enhanced font 'Times,10' size 12cm,9cm
set output '{png_path}'
replot
unset output
"#
    );

    let run_result = run_gnuplot(&script);
    let _ = std::fs::remove_file(dat_path);
    run_result.map(|_| format!("{}, {}", pdf_path, png_path))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_viridis_palette_string_contains_hex_stops() {
        let palette = viridis_palette_string();
        assert!(palette.starts_with("set palette defined ("));
        assert!(palette.contains("0.000 '#440153'"));
        assert!(palette.contains("1.000 '#fde724'"));
    }

    #[test]
    fn test_plot_filename_generation() {
        let (iv_pdf, iv_png) = iv_plot_paths(DisplayMode::Conductance);
        assert_eq!(iv_pdf, "iv_conductance.pdf");
        assert_eq!(iv_png, "iv_conductance.png");

        let (h_pdf, h_png) = heatmap_plot_paths(AppMode::Temperature, DisplayMode::Iets);
        assert_eq!(h_pdf, "temperature_iets.pdf");
        assert_eq!(h_png, "temperature_iets.png");
    }
}
