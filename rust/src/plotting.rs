use std::io::Write;
use std::process::{Command, Stdio};

/// Generate publication-quality I-V plots using gnuplot (PDF + PNG).
pub fn plot_iv(
    csv_path: &str,
    pdf_path: &str,
    png_path: &str,
    n: usize,
    lambda: f64,
    t: f64,
    vmode: f64,
    alpha_l: f64,
    alpha_r: f64,
    eta: f64,
    vg: f64,
    wall_time: f64,
) {
    let script = format!(
        r#"set datafile separator ','
set key left top box opaque
set grid lw 0.5 lc rgb '#cccccc'
set xlabel 'V_{{sd}} (V)'
set ylabel '{{/Italic I}} (A)'
set title 'N={n}, lambda={lambda:.1}, T={t:.1} K, hbar*omega={vmode_mev:.0} meV'
set label 1 'Rust: alpha_L={alpha_l:.2}, alpha_R={alpha_r:.2}, eta={eta:.1}, Vg={vg:.1} V\nt={wall_time:.2} s' at graph 0.98,0.03 right front
set style line 1 lc rgb '#111111' pt 7 ps 0.35 lw 1.0
set style line 2 lc rgb '#0059b3' lw 1.5
set style line 3 lc rgb '#cc1f1f' lw 1.5 dt 2
set terminal pdfcairo enhanced font 'Times,10' size 12cm,9cm
set output '{pdf_path}'
plot '{csv_path}' every ::1 using 1:2 with points ls 1 title 'I_{{total}}', \
     '' every ::1 using 1:3 with lines ls 2 title 'I_{{seq}}', \
     '' every ::1 using 1:4 with lines ls 3 title 'I_{{cot}}'
set terminal pngcairo enhanced font 'Times,10' size 12cm,9cm
set output '{png_path}'
replot
unset output
"#,
        n = n,
        lambda = lambda,
        t = t,
        vmode_mev = vmode * 1e3,
        alpha_l = alpha_l,
        alpha_r = alpha_r,
        eta = eta,
        vg = vg,
        wall_time = wall_time,
        pdf_path = pdf_path,
        csv_path = csv_path,
        png_path = png_path,
    );

    match Command::new("gnuplot")
        .stdin(Stdio::piped())
        .stdout(Stdio::null())
        .stderr(Stdio::inherit())
        .spawn()
    {
        Ok(mut child) => {
            if let Some(stdin) = child.stdin.as_mut() {
                let _ = stdin.write_all(script.as_bytes());
            }
            let status = child.wait();
            if let Ok(s) = status {
                if !s.success() {
                    eprintln!("Warning: gnuplot failed to generate plots.");
                }
            }
        }
        Err(_) => {
            eprintln!("Warning: gnuplot not available, skipping plot generation.");
        }
    }
}
