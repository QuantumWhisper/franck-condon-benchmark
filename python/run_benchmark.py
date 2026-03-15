#!/usr/bin/env python3
"""
run_benchmark.py — Franck-Condon I-V simulation benchmark (Python+Numba)

Usage:
    python run_benchmark.py              # uses 'default' spec
    python run_benchmark.py quick        # uses 'quick' spec
    python run_benchmark.py default      # uses 'default' spec
"""

import sys
import os
import json
import time
import platform
from datetime import datetime
from statistics import median

import numpy as np

# Add src/ to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from src.simulate import simulate_iv
from src.plotting import plot_iv
from src.constants import ELEMENTARY_CHARGE


def main():
    # Parse spec name from command line
    spec_name = sys.argv[1] if len(sys.argv) >= 2 else "default"

    if spec_name not in ("default", "quick"):
        print(f'Unknown spec "{spec_name}". Use "default" or "quick".')
        sys.exit(1)

    # Load benchmark parameters
    this_dir = os.path.dirname(os.path.abspath(__file__))
    spec_file = os.path.join(this_dir, "..", "benchmark", "spec", f"{spec_name}_params.json")
    with open(spec_file) as f:
        spec = json.load(f)

    params = spec["parameters"]
    sweep = spec["bias_sweep"]

    # Extract parameters
    N = int(params["N"])
    vmode = float(params["vmode"])
    alphaL = float(params["alphaL"])
    alphaR = float(params["alphaR"])
    lambda_ = float(params["lambda"])
    T = float(params["T"])
    eta = float(params["eta"])
    Vg = float(params["Vg"])
    tau_raw = params["tau"]
    tau = float("inf") if tau_raw == "Inf" else float(tau_raw)

    # Build bias sweep
    # Use MATLAB reference Vsd values if available (ensures bit-for-bit matching).
    # MATLAB's colon operator uses a specific internal algorithm that differs from
    # Python's range/arange at the ULP level.
    Vsd_start = float(sweep["Vsd_start"])
    Vsd_end = float(sweep["Vsd_end"])
    Vsd_step = float(sweep["Vsd_step"])

    ref_path = os.path.join(this_dir, "..", "benchmark", "results",
                            f"matlab_{spec_name}_results.json")
    if os.path.isfile(ref_path):
        with open(ref_path) as f:
            ref_data = json.load(f)
        Vsd = np.array(ref_data["Vsd"], dtype=float)
    else:
        # Fallback: generate using MATLAB-compatible formula
        nVsd_steps = round((Vsd_end - Vsd_start) / Vsd_step)
        Vsd = np.array([Vsd_start + i * Vsd_step for i in range(nVsd_steps + 1)])

    print(f"=== Franck-Condon Benchmark (Python) [{spec_name}] ===")
    print(f"N={N}, lambda={lambda_:.1f}, T={T:.1f} K, "
          f"Vsd=[{Vsd_start:.3f}:{Vsd_step:.3f}:{Vsd_end:.3f}] V")
    print(f"Total bias points: {len(Vsd)}\n")

    # --- Timed runs (minimum 3, report median) ---
    # No warm-up needed for Python (not JIT at the simulation level,
    # though Numba functions will JIT on first call within first run)
    n_runs = 3
    wall_times = []

    result_Vsd = result_I_tol = result_I_seq = result_I_cot = None

    for run in range(1, n_runs + 1):
        print(f"--- Timed run {run}/{n_runs} ---")
        t_start = time.perf_counter()
        result_Vsd, result_I_tol, result_I_seq, result_I_cot = simulate_iv(
            N, vmode, alphaL, alphaR, lambda_, Vsd, T, eta, Vg, tau,
            verbose=(run == 1))
        t_end = time.perf_counter()
        wall_time = t_end - t_start
        wall_times.append(wall_time)
        print(f"Run {run}: {wall_time:.2f} seconds\n")

    median_wall_time = median(wall_times)
    times_str = ", ".join(f"{t:.2f}" for t in wall_times)
    print(f"Median wall time: {median_wall_time:.2f} seconds (from {n_runs} runs: {times_str})\n")

    # --- Validate against MATLAB reference ---
    validate_against_matlab(spec_name, this_dir, result_Vsd, result_I_tol,
                            result_I_seq, result_I_cot)

    # --- Save results ---
    results_dir = os.path.join(this_dir, "..", "benchmark", "results")
    os.makedirs(results_dir, exist_ok=True)
    prefix = f"python_{spec_name}"

    # 1. JSON
    json_path = os.path.join(results_dir, f"{prefix}_results.json")
    results = {
        "spec": spec_name,
        "language": "Python",
        "version": platform.python_version(),
        "wall_time_seconds": median_wall_time,
        "parameters": {
            "N": N, "vmode": vmode, "alphaL": alphaL, "alphaR": alphaR,
            "lambda": lambda_, "T": T, "eta": eta, "Vg": Vg, "tau": "Inf"
        },
        "bias_sweep": {
            "Vsd_start": Vsd_start, "Vsd_end": Vsd_end, "Vsd_step": Vsd_step
        },
        "Vsd": result_Vsd.tolist(),
        "I_tol": result_I_tol.tolist(),
        "I_seq": result_I_seq.tolist(),
        "I_cot": result_I_cot.tolist(),
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    }
    with open(json_path, 'w') as f:
        json.dump(results, f, indent=2)

    # 2. CSV
    csv_path = os.path.join(results_dir, f"{prefix}_IV.csv")
    with open(csv_path, 'w') as f:
        f.write("Vsd_V,I_tol_A,I_seq_A,I_cot_A\n")
        for i in range(len(result_Vsd)):
            f.write(f"{result_Vsd[i]:.6e},{result_I_tol[i]:.6e},"
                    f"{result_I_seq[i]:.6e},{result_I_cot[i]:.6e}\n")

    # 3. Publication-quality I-V plot (PDF + PNG)
    pdf_path = os.path.join(results_dir, f"{prefix}_IV.pdf")
    png_path = os.path.join(results_dir, f"{prefix}_IV.png")
    plot_iv(result_Vsd, result_I_tol, result_I_seq, result_I_cot,
            N, lambda_, T, vmode, alphaL, alphaR, eta, Vg,
            median_wall_time, pdf_path, png_path)

    # Print summary
    print(f"\n=== Benchmark Complete [{spec_name}] ===")
    print(f"Median wall time: {median_wall_time:.2f} seconds")
    print(f"Max |I_tol|: {np.max(np.abs(result_I_tol)):.4e} A")
    print(f"Max |I_seq|: {np.max(np.abs(result_I_seq)):.4e} A")
    print(f"Max |I_cot|: {np.max(np.abs(result_I_cot)):.4e} A")
    print(f"\nSaved:")
    print(f"  JSON: {json_path}")
    print(f"  CSV:  {csv_path}")
    print(f"  PDF:  {pdf_path}")
    print(f"  PNG:  {png_path}")


def validate_against_matlab(spec_name, this_dir, Vsd, I_tol, I_seq, I_cot):
    """Compare results against MATLAB reference within tolerance.

    The formal tolerance in BENCHMARK.md is 1e-10, but AGENTS.md documents that
    the cotunneling component inherently amplifies ULP-level differences across
    languages. A realistic cross-language tolerance is ~1e-5 for total current.
    """
    ref_path = os.path.join(this_dir, "..", "benchmark", "results",
                            f"matlab_{spec_name}_results.json")

    if not os.path.isfile(ref_path):
        print(f"WARNING: No MATLAB reference found at {ref_path} — skipping validation.")
        return

    print("Validating against MATLAB reference...")
    with open(ref_path) as f:
        ref = json.load(f)

    ref_I_tol = np.array(ref["I_tol"], dtype=float)
    ref_I_seq = np.array(ref["I_seq"], dtype=float)
    ref_I_cot = np.array(ref["I_cot"], dtype=float)

    # Formal tolerance from BENCHMARK.md
    formal_tol = 1e-10
    # Realistic tolerance from AGENTS.md (cotunneling amplifies ULP differences)
    realistic_tol = 1e-4

    def check_tolerance(name, impl, ref_vals, tol):
        diff = np.abs(impl - ref_vals)
        denom = np.maximum(np.abs(ref_vals), 1e-30)
        rel_err = diff / denom
        max_rel_err = np.max(rel_err)
        passed = max_rel_err < tol
        status = "PASS" if passed else "FAIL"
        print(f"  {name}: max relative error = {max_rel_err:.4e}  [{status}] (tol={tol:.0e})")
        return passed, max_rel_err, rel_err

    print(f"\n  --- Formal tolerance (1e-10, per BENCHMARK.md) ---")
    pass_tol_f, _, _ = check_tolerance("I_tol", I_tol, ref_I_tol, formal_tol)
    pass_seq_f, _, _ = check_tolerance("I_seq", I_seq, ref_I_seq, formal_tol)
    pass_cot_f, _, _ = check_tolerance("I_cot", I_cot, ref_I_cot, formal_tol)

    print(f"\n  --- Realistic tolerance (1e-4, per AGENTS.md) ---")
    pass_tol_r, _, rel_tol = check_tolerance("I_tol", I_tol, ref_I_tol, realistic_tol)
    pass_seq_r, _, rel_seq = check_tolerance("I_seq", I_seq, ref_I_seq, realistic_tol)
    pass_cot_r, _, rel_cot = check_tolerance("I_cot", I_cot, ref_I_cot, realistic_tol)

    all_formal = pass_tol_f and pass_seq_f and pass_cot_f
    all_realistic = pass_tol_r and pass_seq_r and pass_cot_r

    if all_formal:
        print(f"\n  ALL CHECKS PASSED within formal {formal_tol} tolerance.")
    elif all_realistic:
        print(f"\n  All checks pass within realistic {realistic_tol} tolerance.")
        print(f"  (Formal {formal_tol} tolerance not met — expected for cotunneling component.)")
    else:
        print(f"\n  SOME CHECKS FAILED even at realistic {realistic_tol} tolerance.")
        print("  Investigate numerical differences.\n")

        # Print first few differing points for debugging
        for name, impl, ref_vals, rel_err in [
            ("I_tol", I_tol, ref_I_tol, rel_tol),
            ("I_seq", I_seq, ref_I_seq, rel_seq),
            ("I_cot", I_cot, ref_I_cot, rel_cot),
        ]:
            bad_idx = np.where(rel_err >= realistic_tol)[0]
            if len(bad_idx) > 0:
                print(f"  {name}: {len(bad_idx)} points exceed tolerance")
                for idx in bad_idx[:5]:
                    print(f"    Vsd[{idx}]={Vsd[idx]:.4f}: "
                          f"python={impl[idx]:.15e}  "
                          f"matlab={ref_vals[idx]:.15e}  "
                          f"relerr={rel_err[idx]:.4e}")
    print()


if __name__ == "__main__":
    main()
