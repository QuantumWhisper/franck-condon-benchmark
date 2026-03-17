use ratatui::{
    layout::{Alignment, Constraint, Layout, Position, Rect},
    style::{Color, Modifier, Style},
    symbols,
    text::Line,
    widgets::{Axis, Block, Chart, Dataset, Gauge, GraphType, LegendPosition, Paragraph},
    Frame,
};

use super::widgets::{Colorbar, Heatmap};
use super::{App, AppMode, DisplayMode};

pub(crate) fn render(frame: &mut Frame, app: &App) {
    let [main_area, help_area] =
        Layout::vertical([Constraint::Fill(1), Constraint::Length(1)]).areas(frame.area());

    let [params_area, plot_area] = Layout::horizontal([
        Constraint::Length(super::PARAM_PANEL_WIDTH),
        Constraint::Fill(1),
    ])
    .areas(main_area);

    render_params(frame, params_area, app);

    match app.mode {
        AppMode::IVCurve => render_iv_chart(frame, plot_area, app),
        AppMode::Stability => render_stability(frame, plot_area, app),
        AppMode::Temperature => render_temperature(frame, plot_area, app),
    }

    render_help(frame, help_area, app);
}

fn render_params(frame: &mut Frame, area: Rect, app: &App) {
    let params = app.visible_params();
    let mode_str = format!("{}{}", app.mode.label(), app.display_mode.suffix());
    let block = Block::bordered().title(Line::from(format!(" {} ", mode_str)).centered());
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
            Style::new().fg(Color::Yellow).add_modifier(Modifier::BOLD)
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
            ("  Press Enter".into(), Style::new().fg(Color::DarkGray))
        };

        let heatmap_mode = app.mode != AppMode::IVCurve;

        let mut lines = vec![text];

        if heatmap_mode && !app.compute_running {
            let mut contrast_parts = Vec::new();
            if let Some(cmin) = app.contrast_min {
                contrast_parts.push(format!("min={:.1}", cmin));
            }
            if let Some(cmax) = app.contrast_max {
                contrast_parts.push(format!("max={:.1}", cmax));
            }
            if (app.gamma - 1.0).abs() > 0.01 {
                contrast_parts.push(format!("γ={:.2}", app.gamma));
            }

            if !contrast_parts.is_empty() {
                lines.push(format!("  Contrast: {}", contrast_parts.join("  ")));
            }
        }

        if app.cursor_mode && heatmap_mode {
            let (x_label, x_unit, x_val, y_val, data_val) = match app.mode {
                AppMode::Stability => {
                    let x = app
                        .stability_vg_vals
                        .get(app.cursor_col)
                        .copied()
                        .unwrap_or(0.0);
                    let y = app
                        .stability_vsd_vals
                        .get(app.cursor_row)
                        .copied()
                        .unwrap_or(0.0);
                    let z = app
                        .stability_grid
                        .get(app.cursor_col)
                        .and_then(|col| col.get(app.cursor_row))
                        .copied()
                        .unwrap_or(0.0);
                    ("Vg", "V", x, y, z)
                }
                AppMode::Temperature => {
                    let x = app
                        .temperature_t_vals
                        .get(app.cursor_col)
                        .copied()
                        .unwrap_or(0.0);
                    let y = app
                        .temperature_vsd_vals
                        .get(app.cursor_row)
                        .copied()
                        .unwrap_or(0.0);
                    let z = app
                        .temperature_grid
                        .get(app.cursor_col)
                        .and_then(|col| col.get(app.cursor_row))
                        .copied()
                        .unwrap_or(0.0);
                    ("T", "K", x, y, z)
                }
                AppMode::IVCurve => ("", "", 0.0, 0.0, 0.0),
            };
            lines.push(format!(
                "  {}={:.3}{}  Vsd={:.3}V",
                x_label, x_val, x_unit, y_val
            ));
            lines.push(format!("  {}={:.3e}", app.display_mode.y_label(), data_val));
        }

        frame.render_widget(Paragraph::new(lines.join("\n")).style(style), status_area);
    }
}

fn render_iv_chart(frame: &mut Frame, area: Rect, app: &App) {
    let dm = app.display_mode;
    let title = format!(
        " {}  N={}  l={:.1}  T={:.1}K ",
        match dm {
            DisplayMode::Current => "I-V Curve",
            DisplayMode::Conductance => "dI/dV",
            DisplayMode::Iets => "IETS (d\u{00b2}I/dV\u{00b2})",
            DisplayMode::NormalizedIets => "Norm. IETS",
        },
        app.n,
        app.lambda,
        app.t_kelvin
    );
    let block = Block::bordered().title(Line::from(title).centered());

    if let Some((ref data, _)) = app.iv_data {
        let vals_tol = dm.transform(&data.vsd, &data.i_tol);
        let vals_seq = dm.transform(&data.vsd, &data.i_seq);
        let vals_cot = dm.transform(&data.vsd, &data.i_cot);

        let data_tol: Vec<(f64, f64)> = data
            .vsd
            .iter()
            .zip(vals_tol.iter())
            .map(|(&x, &y)| (x, y))
            .collect();
        let data_seq: Vec<(f64, f64)> = data
            .vsd
            .iter()
            .zip(vals_seq.iter())
            .map(|(&x, &y)| (x, y))
            .collect();
        let data_cot: Vec<(f64, f64)> = data
            .vsd
            .iter()
            .zip(vals_cot.iter())
            .map(|(&x, &y)| (x, y))
            .collect();

        let y_label = dm.y_label();
        let (n_tol, n_seq, n_cot) = dm.dataset_names();

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
                .name(n_tol)
                .marker(symbols::Marker::Braille)
                .graph_type(GraphType::Line)
                .style(Style::default().fg(Color::Green))
                .data(&data_tol),
            Dataset::default()
                .name(n_seq)
                .marker(symbols::Marker::Braille)
                .graph_type(GraphType::Line)
                .style(Style::default().fg(Color::Yellow))
                .data(&data_seq),
            Dataset::default()
                .name(n_cot)
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
            .title(y_label)
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

fn render_heatmap_crosshair(
    frame: &mut Frame,
    heatmap_area: Rect,
    cursor_col: usize,
    cursor_row: usize,
    n_cols: usize,
    n_rows: usize,
) {
    let dw = heatmap_area.width as usize;
    let dh = heatmap_area.height as usize;
    if n_cols == 0 || n_rows == 0 || dw == 0 || dh == 0 {
        return;
    }

    let clamped_col = cursor_col.min(n_cols - 1);
    let screen_x = if n_cols > 1 && dw > 1 {
        heatmap_area.x + (clamped_col * (dw - 1) / (n_cols - 1)) as u16
    } else {
        heatmap_area.x
    };

    let pixel_h = dh * 2;
    let clamped_row = cursor_row.min(n_rows - 1);
    let target_pixel = if n_rows > 1 && pixel_h > 1 {
        (pixel_h - 1) - clamped_row * (pixel_h - 1) / (n_rows - 1)
    } else {
        0
    };
    let screen_y = heatmap_area.y + (target_pixel / 2) as u16;

    let dim_style = Style::new().fg(Color::DarkGray);

    for dy in 0..heatmap_area.height {
        let y = heatmap_area.y + dy;
        if y != screen_y {
            if let Some(cell) = frame.buffer_mut().cell_mut(Position::new(screen_x, y)) {
                cell.set_symbol("│").set_style(dim_style);
            }
        }
    }

    for dx in 0..heatmap_area.width {
        let x = heatmap_area.x + dx;
        if x != screen_x {
            if let Some(cell) = frame.buffer_mut().cell_mut(Position::new(x, screen_y)) {
                cell.set_symbol("─").set_style(dim_style);
            }
        }
    }

    if let Some(cell) = frame
        .buffer_mut()
        .cell_mut(Position::new(screen_x, screen_y))
    {
        cell.set_symbol("┼")
            .set_style(Style::new().fg(Color::White).add_modifier(Modifier::BOLD));
    }
}

fn render_line_profile(frame: &mut Frame, area: Rect, app: &App, display_grid: &[Vec<f64>]) {
    let n_cols = display_grid.len();
    let n_rows = display_grid.iter().map(|r| r.len()).max().unwrap_or(0);
    if n_cols == 0 || n_rows == 0 {
        return;
    }

    if app.line_profile_horizontal {
        let col = app.cursor_col.min(n_cols - 1);
        let row_data = &display_grid[col];
        if row_data.is_empty() {
            return;
        }

        let (x_vals, x_label) = match app.mode {
            AppMode::Stability => (&app.stability_vsd_vals, "Vsd (V)"),
            AppMode::Temperature => (&app.temperature_vsd_vals, "Vsd (V)"),
            AppMode::IVCurve => return,
        };

        let (col_label, col_val) = match app.mode {
            AppMode::Stability => ("Vg", app.stability_vg_vals.get(col).copied().unwrap_or(0.0)),
            AppMode::Temperature => ("T", app.temperature_t_vals.get(col).copied().unwrap_or(0.0)),
            AppMode::IVCurve => return,
        };

        let data_points: Vec<(f64, f64)> = x_vals
            .iter()
            .zip(row_data.iter())
            .map(|(&x, &y)| (x, y))
            .collect();

        render_profile_chart(
            frame,
            area,
            &data_points,
            x_label,
            &format!(
                "{} @ {}={:.3}",
                app.display_mode.y_label(),
                col_label,
                col_val
            ),
        );
    } else {
        let row = app.cursor_row.min(n_rows - 1);
        let cut_data: Vec<f64> = display_grid
            .iter()
            .map(|col| col.get(row).copied().unwrap_or(0.0))
            .collect();

        let (x_vals, x_label, row_val) = match app.mode {
            AppMode::Stability => (
                &app.stability_vg_vals,
                "Vg (V)",
                app.stability_vsd_vals.get(row).copied().unwrap_or(0.0),
            ),
            AppMode::Temperature => (
                &app.temperature_t_vals,
                "T (K)",
                app.temperature_vsd_vals.get(row).copied().unwrap_or(0.0),
            ),
            AppMode::IVCurve => return,
        };

        let data_points: Vec<(f64, f64)> = x_vals
            .iter()
            .zip(cut_data.iter())
            .map(|(&x, &y)| (x, y))
            .collect();

        render_profile_chart(
            frame,
            area,
            &data_points,
            x_label,
            &format!("{} @ Vsd={:.3}V", app.display_mode.y_label(), row_val),
        );
    }
}

fn render_profile_chart(
    frame: &mut Frame,
    area: Rect,
    data: &[(f64, f64)],
    x_label: &str,
    title: &str,
) {
    if data.is_empty() {
        return;
    }

    let block = Block::bordered().title(Line::from(format!(" {} ", title)).centered());

    let x_min = data.first().map(|d| d.0).unwrap_or(0.0);
    let x_max = data.last().map(|d| d.0).unwrap_or(1.0);
    let y_min = data.iter().map(|d| d.1).fold(f64::INFINITY, f64::min);
    let y_max = data.iter().map(|d| d.1).fold(f64::NEG_INFINITY, f64::max);
    let pad = (y_max - y_min).abs().max(1e-30) * 0.08;
    let y_lo = y_min - pad;
    let y_hi = y_max + pad;

    let datasets = vec![Dataset::default()
        .name("cut")
        .marker(symbols::Marker::Braille)
        .graph_type(GraphType::Line)
        .style(Style::default().fg(Color::Green))
        .data(data)];

    let x_axis = Axis::default()
        .title(x_label)
        .style(Style::default().fg(Color::Gray))
        .bounds([x_min, x_max])
        .labels(vec![format!("{:.2}", x_min), format!("{:.2}", x_max)]);

    let y_axis = Axis::default()
        .style(Style::default().fg(Color::Gray))
        .bounds([y_lo, y_hi])
        .labels(vec![format!("{:.1e}", y_lo), format!("{:.1e}", y_hi)]);

    let chart = Chart::new(datasets)
        .block(block)
        .x_axis(x_axis)
        .y_axis(y_axis);

    frame.render_widget(chart, area);
}

fn render_stability(frame: &mut Frame, area: Rect, app: &App) {
    let title = format!(
        " Stability{}  N={}  l={:.1}  T={:.1}K ",
        app.display_mode.suffix(),
        app.n,
        app.lambda,
        app.t_kelvin
    );

    let prog_h = if app.compute_running || app.stability_progress.1 > 0 {
        3u16
    } else {
        0
    };
    let [map_area, progress_area] =
        Layout::vertical([Constraint::Fill(1), Constraint::Length(prog_h)]).areas(area);

    let (heatmap_outer, profile_area) = if app.cursor_mode && app.show_line_profile {
        let [top, bottom] =
            Layout::vertical([Constraint::Percentage(65), Constraint::Percentage(35)])
                .areas(map_area);
        (top, Some(bottom))
    } else {
        (map_area, None)
    };

    let block = Block::bordered().title(Line::from(title).centered());
    let has_data = app.stability_grid.iter().any(|r| !r.is_empty());

    if has_data {
        let inner = block.inner(heatmap_outer);
        frame.render_widget(block, heatmap_outer);

        let grid_buf: Vec<Vec<f64>>;
        let display_grid: &[Vec<f64>] = if app.display_mode == DisplayMode::Current {
            &app.stability_grid
        } else {
            grid_buf = app
                .display_mode
                .transform_grid(&app.stability_grid, &app.stability_vsd_vals);
            &grid_buf
        };

        let mut log_min = f64::INFINITY;
        let mut log_max = f64::NEG_INFINITY;
        for row in display_grid {
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
        let effective_min = app.contrast_min.unwrap_or(log_min);
        let effective_max = app.contrast_max.unwrap_or(log_max);

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
                data: display_grid,
                log_min: effective_min,
                log_max: effective_max,
                gamma: app.gamma,
            },
            heatmap_area,
        );

        if app.cursor_mode {
            let n_cols = display_grid.len();
            let n_rows = display_grid.iter().map(|r| r.len()).max().unwrap_or(0);
            render_heatmap_crosshair(
                frame,
                heatmap_area,
                app.cursor_col,
                app.cursor_row,
                n_cols,
                n_rows,
            );
        }

        frame.render_widget(Colorbar { gamma: app.gamma }, colorbar_area);

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
                Paragraph::new(format!("{:.1}", effective_max)).style(label_style),
                cb_top,
            );
            frame.render_widget(
                Paragraph::new(format!("{:.1}", effective_min)).style(label_style),
                cb_bot,
            );
            if cblabel_area.height >= 4 {
                let cb_title = Rect::new(cblabel_area.x, cblabel_area.y + 1, cblabel_area.width, 1);
                frame.render_widget(
                    Paragraph::new(app.display_mode.colorbar_label())
                        .style(Style::new().fg(Color::DarkGray)),
                    cb_title,
                );
            }
        }

        if let Some(prof_area) = profile_area {
            render_line_profile(frame, prof_area, app, display_grid);
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
            heatmap_outer,
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

fn render_temperature(frame: &mut Frame, area: Rect, app: &App) {
    let title = format!(
        " Temperature{}  N={}  l={:.1}  Vg={:.3}V ",
        app.display_mode.suffix(),
        app.n,
        app.lambda,
        app.vg
    );

    let prog_h = if app.compute_running || app.temperature_progress.1 > 0 {
        3u16
    } else {
        0
    };
    let [map_area, progress_area] =
        Layout::vertical([Constraint::Fill(1), Constraint::Length(prog_h)]).areas(area);

    let (heatmap_outer, profile_area) = if app.cursor_mode && app.show_line_profile {
        let [top, bottom] =
            Layout::vertical([Constraint::Percentage(65), Constraint::Percentage(35)])
                .areas(map_area);
        (top, Some(bottom))
    } else {
        (map_area, None)
    };

    let block = Block::bordered().title(Line::from(title).centered());
    let has_data = app.temperature_grid.iter().any(|r| !r.is_empty());

    if has_data {
        let inner = block.inner(heatmap_outer);
        frame.render_widget(block, heatmap_outer);

        let grid_buf: Vec<Vec<f64>>;
        let display_grid: &[Vec<f64>] = if app.display_mode == DisplayMode::Current {
            &app.temperature_grid
        } else {
            grid_buf = app
                .display_mode
                .transform_grid(&app.temperature_grid, &app.temperature_vsd_vals);
            &grid_buf
        };

        let mut log_min = f64::INFINITY;
        let mut log_max = f64::NEG_INFINITY;
        for row in display_grid {
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
        let effective_min = app.contrast_min.unwrap_or(log_min);
        let effective_max = app.contrast_max.unwrap_or(log_max);

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
                data: display_grid,
                log_min: effective_min,
                log_max: effective_max,
                gamma: app.gamma,
            },
            heatmap_area,
        );

        if app.cursor_mode {
            let n_cols = display_grid.len();
            let n_rows = display_grid.iter().map(|r| r.len()).max().unwrap_or(0);
            render_heatmap_crosshair(
                frame,
                heatmap_area,
                app.cursor_col,
                app.cursor_row,
                n_cols,
                n_rows,
            );
        }

        frame.render_widget(Colorbar { gamma: app.gamma }, colorbar_area);

        let vsd_max_val = app
            .temperature_vsd_vals
            .last()
            .copied()
            .unwrap_or(app.temp_vsd_max);
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

        let t_lo = app.temperature_t_vals.first().copied().unwrap_or(app.t_min);
        let t_hi = app.temperature_t_vals.last().copied().unwrap_or(app.t_max);
        let t_mid = (t_lo + t_hi) / 2.0;
        let xl = Rect::new(heatmap_area.x, xlabel_row.y, heatmap_area.width, 1);

        if xl.width >= 20 {
            let s_left = format!("{:.1}K", t_lo);
            let s_mid = format!("{:.1}K", t_mid);
            let s_right = format!("{:.1}K", t_hi);

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
                Paragraph::new(format!("{:.1}", effective_max)).style(label_style),
                cb_top,
            );
            frame.render_widget(
                Paragraph::new(format!("{:.1}", effective_min)).style(label_style),
                cb_bot,
            );
            if cblabel_area.height >= 4 {
                let cb_title = Rect::new(cblabel_area.x, cblabel_area.y + 1, cblabel_area.width, 1);
                frame.render_widget(
                    Paragraph::new(app.display_mode.colorbar_label())
                        .style(Style::new().fg(Color::DarkGray)),
                    cb_title,
                );
            }
        }

        if let Some(prof_area) = profile_area {
            render_line_profile(frame, prof_area, app, display_grid);
        }
    } else {
        let msg = if app.compute_running {
            "Computing temperature diagram..."
        } else {
            "Press Enter to compute temperature diagram"
        };
        frame.render_widget(
            Paragraph::new(msg)
                .block(block)
                .style(Style::default().fg(Color::DarkGray)),
            heatmap_outer,
        );
    }

    if prog_h > 0 {
        let (done, total) = app.temperature_progress;
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
    } else if app.cursor_mode && app.mode != AppMode::IVCurve {
        " Arrows Move | h Profile | v H/V cut | d Display | c Exit cursor | e Csv | p Plot | s Save | q Quit "
    } else {
        match app.mode {
            AppMode::IVCurve => {
                " Up/Dn Select | Lt/Rt Adj (Shift=fine) auto | d I/G/IETS/nIETS | Tab Mode | e Csv | p Plot | s Save | q Quit "
            }
            _ => {
                " Up/Dn Select | Lt/Rt Adj | Enter Run | d Display | [/] Range | g/G Gamma | a Auto | c Cursor | Tab | e Csv | p Plot | s Save | q Quit "
            }
        }
    };
    frame.render_widget(
        Paragraph::new(help).style(Style::default().fg(Color::DarkGray)),
        area,
    );
}
