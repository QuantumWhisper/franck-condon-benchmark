use ratatui::{
    buffer::Buffer,
    layout::Rect,
    style::{Color, Style},
    widgets::Widget,
};

pub(crate) const VIRIDIS: [(f64, f64, f64); 9] = [
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

pub(crate) fn viridis(t: f64) -> Color {
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

pub(crate) struct Heatmap<'a> {
    pub(crate) data: &'a [Vec<f64>],
    pub(crate) log_min: f64,
    pub(crate) log_max: f64,
    pub(crate) gamma: f64,
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

                let top_t = ((top_val - self.log_min) / log_range)
                    .clamp(0.0, 1.0)
                    .powf(self.gamma);
                let bot_t = ((bot_val - self.log_min) / log_range)
                    .clamp(0.0, 1.0)
                    .powf(self.gamma);

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

pub(crate) struct Colorbar {
    pub(crate) gamma: f64,
}

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
                ((pixel_h - 1 - top_pixel) as f64 / (pixel_h - 1) as f64).powf(self.gamma)
            } else {
                0.5
            };

            let bot_pixel = dy * 2 + 1;
            let bot_t = if bot_pixel < pixel_h && pixel_h > 1 {
                ((pixel_h - 1 - bot_pixel) as f64 / (pixel_h - 1) as f64).powf(self.gamma)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_viridis_start() {
        if let Color::Rgb(r, g, b) = viridis(0.0) {
            assert_eq!(r, 68);
            assert_eq!(g, 1);
            assert_eq!(b, 83);
        } else {
            panic!("Expected Rgb");
        }
    }

    #[test]
    fn test_viridis_end() {
        if let Color::Rgb(r, g, b) = viridis(1.0) {
            assert_eq!(r, 253);
            assert_eq!(g, 231);
            assert_eq!(b, 36);
        } else {
            panic!("Expected Rgb");
        }
    }

    #[test]
    fn test_viridis_clamp_negative() {
        assert_eq!(viridis(-1.0), viridis(0.0));
    }

    #[test]
    fn test_viridis_clamp_above() {
        assert_eq!(viridis(2.0), viridis(1.0));
    }

    #[test]
    fn test_gamma_identity() {
        let t: f64 = 0.5;
        assert!((t.powf(1.0) - t).abs() < 1e-15);
    }

    #[test]
    fn test_gamma_brightens() {
        let t: f64 = 0.25;
        assert!(t.powf(0.5) > t);
    }
}
