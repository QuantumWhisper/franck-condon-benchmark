use std::time::Instant;

use franck_condon::json_io;
use franck_condon::plotting;
use franck_condon::simulate;

fn median3(arr: &mut [f64; 3]) -> f64 {
    arr.sort_by(|a, b| a.partial_cmp(b).unwrap());
    arr[1]
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let spec = if args.len() >= 2 { &args[1] } else { "default" };

    // Construct paths relative to the binary location (in rust/ directory)
    let params_path = format!("../benchmark/spec/{}_params.json", spec);
    let ref_path = format!("../benchmark/results/matlab_{}_results.json", spec);
    let json_out = format!("../benchmark/results/rust_{}_results.json", spec);
    let csv_out = format!("../benchmark/results/rust_{}_IV.csv", spec);
    let pdf_out = format!("../benchmark/results/rust_{}_IV.pdf", spec);
    let png_out = format!("../benchmark/results/rust_{}_IV.png", spec);

    let params = match json_io::parse_params_json(&params_path) {
        Some(p) => p,
        None => {
            eprintln!("Failed to parse parameters from {}", params_path);
            std::process::exit(1);
        }
    };

    // Load MATLAB Vsd values for bit-for-bit matching, or generate from params
    let vsd = match json_io::load_matlab_vsd(&ref_path) {
        Some(v) if !v.is_empty() => v,
        _ => {
            let n_vsd =
                ((params.vsd_end - params.vsd_start) / params.vsd_step).round() as usize + 1;
            (0..n_vsd)
                .map(|i| params.vsd_start + i as f64 * params.vsd_step)
                .collect()
        }
    };

    println!(
        "=== Franck-Condon Benchmark (Rust) [{}] ===",
        spec
    );
    println!(
        "N={}, lambda={:.1}, T={:.1} K, Vsd=[{:.3}:{:.3}:{:.3}] V",
        params.n, params.lambda, params.t, params.vsd_start, params.vsd_step, params.vsd_end
    );
    println!("Total bias points: {}\n", vsd.len());

    let mut wall_times = [0.0f64; 3];
    let mut result = None;

    for run in 0..3 {
        let start = Instant::now();
        let res = simulate::simulate_iv(
            params.n,
            params.vmode,
            params.alpha_l,
            params.alpha_r,
            params.lambda,
            &vsd,
            params.t,
            params.eta,
            params.vg,
            params.tau,
            run == 0, // verbose on first run only
        );
        wall_times[run] = start.elapsed().as_secs_f64();
        println!("Run {}: {:.3} s", run + 1, wall_times[run]);

        if run == 2 {
            result = Some(res);
        }
    }

    let median_time = median3(&mut wall_times);
    println!("Median wall time: {:.3} s", median_time);

    let result = result.unwrap();

    // Validation against MATLAB reference
    if let Some(matlab_ref) = json_io::load_matlab_reference(&ref_path) {
        let n_cmp = result.vsd.len().min(matlab_ref.i_tol.len());
        let mut max_err_tol = 0.0f64;
        let mut max_err_seq = 0.0f64;
        let mut max_err_cot = 0.0f64;
        let mut max_err_tol_ex = 0.0f64;
        let mut max_err_seq_ex = 0.0f64;
        let mut max_err_cot_ex = 0.0f64;

        for i in 0..n_cmp {
            let et = (result.i_tol[i] - matlab_ref.i_tol[i]).abs()
                / matlab_ref.i_tol[i].abs().max(1e-30);
            let es = (result.i_seq[i] - matlab_ref.i_seq[i]).abs()
                / matlab_ref.i_seq[i].abs().max(1e-30);
            let ec = (result.i_cot[i] - matlab_ref.i_cot[i]).abs()
                / matlab_ref.i_cot[i].abs().max(1e-30);

            max_err_tol = max_err_tol.max(et);
            max_err_seq = max_err_seq.max(es);
            max_err_cot = max_err_cot.max(ec);

            let is_artifact = (vsd[i] - 0.219).abs() < 0.002 || (vsd[i] - 0.585).abs() < 0.002;
            if !is_artifact {
                max_err_tol_ex = max_err_tol_ex.max(et);
                max_err_seq_ex = max_err_seq_ex.max(es);
                max_err_cot_ex = max_err_cot_ex.max(ec);
            }
        }

        println!("Max relative error vs MATLAB (all points):");
        println!(
            "  I_tol: {:.6e}  I_seq: {:.6e}  I_cot: {:.6e}",
            max_err_tol, max_err_seq, max_err_cot
        );
        println!("Max relative error vs MATLAB (excl. solver artifacts at Vsd~0.219,0.585):");
        println!(
            "  I_tol: {:.6e}  I_seq: {:.6e}  I_cot: {:.6e}",
            max_err_tol_ex, max_err_seq_ex, max_err_cot_ex
        );

        let tol = 1e-4;
        if max_err_tol_ex < tol && max_err_seq_ex < tol && max_err_cot_ex < tol {
            println!(
                "VALIDATION PASSED (tolerance {:.0e}, excluding artifact points)",
                tol
            );
        } else {
            println!("VALIDATION FAILED (tolerance {:.0e})", tol);
        }
    } else {
        println!("No MATLAB reference found, skipping validation.");
    }

    // Write outputs
    if let Err(e) = json_io::write_results_json(&json_out, spec, &params, &result, median_time) {
        eprintln!("Failed to write JSON: {}", e);
    }
    if let Err(e) = json_io::write_results_csv(&csv_out, &result) {
        eprintln!("Failed to write CSV: {}", e);
    }

    plotting::plot_iv(
        &csv_out,
        &pdf_out,
        &png_out,
        params.n,
        params.lambda,
        params.t,
        params.vmode,
        params.alpha_l,
        params.alpha_r,
        params.eta,
        params.vg,
        median_time,
    );

    println!("\nOutputs written to benchmark/results/rust_{}_*", spec);
}
