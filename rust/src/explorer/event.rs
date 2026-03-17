use std::time::{Duration, Instant};

use ratatui::crossterm::event::{self, KeyCode, KeyModifiers};

use super::export;
use super::{App, AppMode};

pub(crate) fn handle_key(app: &mut App, key: event::KeyEvent) {
    match key.code {
        KeyCode::Char('q') | KeyCode::Char('Q') => {
            app.should_quit = true;
        }
        KeyCode::Esc => {
            if app.cursor_mode {
                app.cursor_mode = false;
                app.show_line_profile = false;
                app.status_msg = Some(("Cursor OFF".into(), false));
            } else if app.compute_running {
                app.cancel_compute();
            } else {
                app.should_quit = true;
            }
        }
        KeyCode::Tab => {
            if !app.compute_running {
                app.cursor_mode = false;
                app.show_line_profile = false;
                app.mode = match app.mode {
                    AppMode::IVCurve => AppMode::Stability,
                    AppMode::Stability => AppMode::Temperature,
                    AppMode::Temperature => AppMode::IVCurve,
                };
                app.selected_param = 0;
                app.iv_recompute_at = None;
            }
        }
        KeyCode::Up => {
            if app.cursor_mode && app.mode != AppMode::IVCurve {
                let max_row = match app.mode {
                    AppMode::Stability => app.stability_vsd_vals.len().saturating_sub(1),
                    AppMode::Temperature => app.temperature_vsd_vals.len().saturating_sub(1),
                    AppMode::IVCurve => 0,
                };
                if app.cursor_row < max_row {
                    app.cursor_row += 1;
                }
            } else if app.selected_param > 0 {
                app.selected_param -= 1;
            }
        }
        KeyCode::Down => {
            if app.cursor_mode && app.mode != AppMode::IVCurve {
                if app.cursor_row > 0 {
                    app.cursor_row -= 1;
                }
            } else {
                let max = app.visible_params().len().saturating_sub(1);
                if app.selected_param < max {
                    app.selected_param += 1;
                }
            }
        }
        KeyCode::Right | KeyCode::Left => {
            if (app.mode == AppMode::Stability || app.mode == AppMode::Temperature)
                && app.compute_running
            {
                return;
            }
            if app.cursor_mode && app.mode != AppMode::IVCurve {
                if key.code == KeyCode::Right {
                    let max_col = match app.mode {
                        AppMode::Stability => app.stability_vg_vals.len().saturating_sub(1),
                        AppMode::Temperature => app.temperature_t_vals.len().saturating_sub(1),
                        AppMode::IVCurve => 0,
                    };
                    if app.cursor_col < max_col {
                        app.cursor_col += 1;
                    }
                } else if app.cursor_col > 0 {
                    app.cursor_col -= 1;
                }
            } else {
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
        }
        KeyCode::Enter => {
            app.iv_recompute_at = None;
            if app.compute_running {
                app.cancel_compute();
            } else {
                match app.mode {
                    AppMode::IVCurve => app.start_iv(),
                    AppMode::Stability => app.start_stability(),
                    AppMode::Temperature => app.start_temperature(),
                }
            }
        }
        KeyCode::Char('[') => {
            if app.mode != AppMode::IVCurve && !app.compute_running {
                let current = app.contrast_max.unwrap_or(-5.0);
                app.contrast_max = Some(current - 1.0);
                app.status_msg = Some((format!("max={:.1}", app.contrast_max.unwrap()), false));
            }
        }
        KeyCode::Char(']') => {
            if app.mode != AppMode::IVCurve && !app.compute_running {
                let current = app.contrast_max.unwrap_or(-5.0);
                app.contrast_max = Some(current + 1.0);
                app.status_msg = Some((format!("max={:.1}", app.contrast_max.unwrap()), false));
            }
        }
        KeyCode::Char('{') => {
            if app.mode != AppMode::IVCurve && !app.compute_running {
                let current = app.contrast_min.unwrap_or(-30.0);
                app.contrast_min = Some(current + 1.0);
                app.status_msg = Some((format!("min={:.1}", app.contrast_min.unwrap()), false));
            }
        }
        KeyCode::Char('}') => {
            if app.mode != AppMode::IVCurve && !app.compute_running {
                let current = app.contrast_min.unwrap_or(-30.0);
                app.contrast_min = Some(current - 1.0);
                app.status_msg = Some((format!("min={:.1}", app.contrast_min.unwrap()), false));
            }
        }
        KeyCode::Char('g') => {
            if app.mode != AppMode::IVCurve && !app.compute_running {
                app.gamma = (app.gamma - 0.1).max(0.1);
                app.status_msg = Some((format!("γ={:.2}", app.gamma), false));
            }
        }
        KeyCode::Char('G') => {
            if app.mode != AppMode::IVCurve && !app.compute_running {
                app.gamma = (app.gamma + 0.1).min(5.0);
                app.status_msg = Some((format!("γ={:.2}", app.gamma), false));
            }
        }
        KeyCode::Char('a') | KeyCode::Char('A') => {
            if app.mode != AppMode::IVCurve && !app.compute_running {
                app.contrast_min = None;
                app.contrast_max = None;
                app.gamma = 1.0;
                app.status_msg = Some(("Auto contrast".into(), false));
            }
        }
        KeyCode::Char('c') | KeyCode::Char('C') => {
            if app.mode != AppMode::IVCurve {
                app.cursor_mode = !app.cursor_mode;
                if !app.cursor_mode {
                    app.show_line_profile = false;
                }
                app.status_msg = Some((
                    if app.cursor_mode {
                        "Cursor ON".into()
                    } else {
                        "Cursor OFF".into()
                    },
                    false,
                ));
            }
        }
        KeyCode::Char('h') | KeyCode::Char('H') => {
            if app.cursor_mode {
                app.show_line_profile = !app.show_line_profile;
            }
        }
        KeyCode::Char('v') | KeyCode::Char('V') => {
            if app.cursor_mode && app.show_line_profile {
                app.line_profile_horizontal = !app.line_profile_horizontal;
            }
        }
        KeyCode::Char('p') | KeyCode::Char('P') => {
            if !app.compute_running {
                let result = match app.mode {
                    AppMode::IVCurve => export::export_iv_gnuplot(app),
                    AppMode::Stability | AppMode::Temperature => {
                        export::export_heatmap_gnuplot(app)
                    }
                };
                match result {
                    Ok(paths) => app.status_msg = Some((format!("Plot: {}", paths), false)),
                    Err(e) => app.status_msg = Some((format!("Plot failed: {}", e), true)),
                }
            }
        }
        KeyCode::Char('d') | KeyCode::Char('D') => {
            app.display_mode = app.display_mode.next();
        }
        KeyCode::Char('e') | KeyCode::Char('E') => {
            if !app.compute_running {
                let result = match app.mode {
                    AppMode::Stability => export::export_stability(app),
                    AppMode::IVCurve => export::export_iv(app),
                    AppMode::Temperature => export::export_temperature(app),
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
