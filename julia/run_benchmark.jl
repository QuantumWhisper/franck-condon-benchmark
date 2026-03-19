#!/usr/bin/env julia
#
# run_benchmark.jl — Franck-Condon I-V simulation benchmark (Julia)
#
# Usage:
#   julia --threads=auto run_benchmark.jl              # uses 'default' spec
#   julia --threads=auto run_benchmark.jl quick        # uses 'quick' spec
#   julia --threads=auto run_benchmark.jl default      # uses 'default' spec
#
# Threading: Use --threads=auto (or -t auto) for parallel bias-point computation.
# Without it, Julia runs single-threaded and misses the ~6-8x parallelism speedup.
#

using Pkg
Pkg.activate(@__DIR__)

using Printf
using Dates
using JSON3
using Statistics

# Include the module
include(joinpath(@__DIR__, "src", "FranckCondon.jl"))
using .FranckCondon

# Include plotting separately (not part of core physics module)
include(joinpath(@__DIR__, "src", "plotting.jl"))

function main()
    # Parse spec name from command line
    spec_name = length(ARGS) >= 1 ? ARGS[1] : "default"

    if spec_name ∉ ("default", "quick")
        error("Unknown spec \"$spec_name\". Use \"default\" or \"quick\".")
    end

    # Load benchmark parameters
    spec_file = joinpath(@__DIR__, "..", "benchmark", "spec", "$(spec_name)_params.json")
    spec = JSON3.read(read(spec_file, String))
    params = spec.parameters
    sweep = spec.bias_sweep

    # Extract parameters
    N = Int(params.N)
    vmode = Float64(params.vmode)
    alphaL = Float64(params.alphaL)
    alphaR = Float64(params.alphaR)
    lambda = Float64(params.lambda)
    T = Float64(params.T)
    eta = Float64(params.eta)
    Vg = Float64(params.Vg)
    tau_raw = params.tau
    tau = (tau_raw isa String && tau_raw == "Inf") ? Inf : Float64(tau_raw)

    # Build bias sweep
    # Use MATLAB reference Vsd values if available (ensures bit-for-bit matching).
    # MATLAB's colon operator uses a specific internal algorithm that differs from
    # Julia's range at the ULP level, making independent generation unreliable.
    Vsd_start = Float64(sweep.Vsd_start)
    Vsd_end = Float64(sweep.Vsd_end)
    Vsd_step = Float64(sweep.Vsd_step)
    ref_path = joinpath(@__DIR__, "..", "benchmark", "results", "matlab_$(spec_name)_results.json")
    if isfile(ref_path)
        ref_data = JSON3.read(read(ref_path, String))
        Vsd = Float64.(ref_data.Vsd)
    else
        # Fallback: generate using MATLAB-compatible formula
        nVsd_steps = round(Int, (Vsd_end - Vsd_start) / Vsd_step)
        Vsd = [Vsd_start + i * Vsd_step for i in 0:nVsd_steps]
    end

    println("=== Franck-Condon Benchmark (Julia) [$spec_name] ===")
    @printf("N=%d, lambda=%.1f, T=%.1f K, Vsd=[%.3f:%.3f:%.3f] V\n",
            N, lambda, T, Vsd_start, Vsd_step, Vsd_end)
    @printf("Total bias points: %d\n", length(Vsd))
    @printf("Julia threads: %d\n", Threads.nthreads())
    if Threads.nthreads() == 1
        println("WARNING: Running single-threaded. Use julia --threads=auto for parallel execution.")
    end
    println()

    # --- Warm-up run (for JIT compilation) ---
    println("Warm-up run (JIT compilation)...")
    simulate_iv(N, vmode, alphaL, alphaR, lambda, Vsd, T, eta, Vg, tau; verbose=false)
    println("Warm-up complete.\n")

    # --- Timed runs (minimum 3, report median) ---
    n_runs = 3
    wall_times = Float64[]

    local result_Vsd, result_I_tol, result_I_seq, result_I_cot

    for run in 1:n_runs
        println("--- Timed run $run/$n_runs ---")
        t_start = time()
        result_Vsd, result_I_tol, result_I_seq, result_I_cot =
            simulate_iv(N, vmode, alphaL, alphaR, lambda, Vsd, T, eta, Vg, tau; verbose=(run==1))
        t_end = time()
        wall_time = t_end - t_start
        push!(wall_times, wall_time)
        @printf("Run %d: %.2f seconds\n\n", run, wall_time)
    end

    median_wall_time = median(wall_times)
    @printf("Median wall time: %.2f seconds (from %d runs: %s)\n\n",
            median_wall_time, n_runs,
            join([@sprintf("%.2f", t) for t in wall_times], ", "))

    # --- Validate against MATLAB reference ---
    validate_against_matlab(spec_name, result_Vsd, result_I_tol, result_I_seq, result_I_cot)

    # --- Save results ---
    results_dir = joinpath(@__DIR__, "..", "benchmark", "results")
    mkpath(results_dir)
    prefix = "julia_$(spec_name)"

    # 1. JSON
    json_path = joinpath(results_dir, "$(prefix)_results.json")
    results = Dict(
        "spec" => spec_name,
        "language" => "Julia",
        "version" => string(VERSION),
        "wall_time_seconds" => median_wall_time,
        "parameters" => Dict(
            "N" => N, "vmode" => vmode, "alphaL" => alphaL, "alphaR" => alphaR,
            "lambda" => lambda, "T" => T, "eta" => eta, "Vg" => Vg, "tau" => "Inf"
        ),
        "bias_sweep" => Dict(
            "Vsd_start" => Vsd_start, "Vsd_end" => Vsd_end, "Vsd_step" => Vsd_step
        ),
        "Vsd" => result_Vsd,
        "I_tol" => result_I_tol,
        "I_seq" => result_I_seq,
        "I_cot" => result_I_cot,
        "timestamp" => Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    )
    open(json_path, "w") do io
        JSON3.pretty(io, results)
    end

    # 2. CSV
    csv_path = joinpath(results_dir, "$(prefix)_IV.csv")
    open(csv_path, "w") do io
        println(io, "Vsd_V,I_tol_A,I_seq_A,I_cot_A")
        for i in eachindex(result_Vsd)
            @printf(io, "%.6e,%.6e,%.6e,%.6e\n",
                    result_Vsd[i], result_I_tol[i], result_I_seq[i], result_I_cot[i])
        end
    end

    # 3. Publication-quality I-V plot (PDF + PNG)
    pdf_path = joinpath(results_dir, "$(prefix)_IV.pdf")
    png_path = joinpath(results_dir, "$(prefix)_IV.png")
    plot_iv(result_Vsd, result_I_tol, result_I_seq, result_I_cot,
            N, lambda, T, vmode, alphaL, alphaR, eta, Vg,
            median_wall_time, pdf_path, png_path)

    # Print summary
    println("\n=== Benchmark Complete [$spec_name] ===")
    @printf("Median wall time: %.2f seconds\n", median_wall_time)
    @printf("Max |I_tol|: %.4e A\n", maximum(abs.(result_I_tol)))
    @printf("Max |I_seq|: %.4e A\n", maximum(abs.(result_I_seq)))
    @printf("Max |I_cot|: %.4e A\n", maximum(abs.(result_I_cot)))
    println("\nSaved:")
    println("  JSON: $json_path")
    println("  CSV:  $csv_path")
    println("  PDF:  $pdf_path")
    println("  PNG:  $png_path")
end

"""
    validate_against_matlab(spec_name, Vsd, I_tol, I_seq, I_cot)

Compare results against MATLAB reference within 1e-10 relative tolerance.
"""
function validate_against_matlab(spec_name, Vsd, I_tol, I_seq, I_cot)
    ref_path = joinpath(@__DIR__, "..", "benchmark", "results", "matlab_$(spec_name)_results.json")

    if !isfile(ref_path)
        println("WARNING: No MATLAB reference found at $ref_path — skipping validation.")
        return
    end

    println("Validating against MATLAB reference...")
    ref = JSON3.read(read(ref_path, String))

    ref_I_tol = Float64.(ref.I_tol)
    ref_I_seq = Float64.(ref.I_seq)
    ref_I_cot = Float64.(ref.I_cot)

    tol = 1e-10

    # Comparison formula from BENCHMARK.md:
    # max(abs(I_impl - I_ref) ./ max(abs(I_ref), 1e-30)) < 1e-10
    function check_tolerance(name, impl, ref_vals)
        diff = abs.(impl .- ref_vals)
        denom = max.(abs.(ref_vals), 1e-30)
        rel_err = diff ./ denom
        max_rel_err = maximum(rel_err)
        passed = max_rel_err < tol
        status = passed ? "PASS" : "FAIL"
        @printf("  %s: max relative error = %.4e  [%s]\n", name, max_rel_err, status)
        return passed
    end

    pass_tol = check_tolerance("I_tol", I_tol, ref_I_tol)
    pass_seq = check_tolerance("I_seq", I_seq, ref_I_seq)
    pass_cot = check_tolerance("I_cot", I_cot, ref_I_cot)

    all_pass = pass_tol && pass_seq && pass_cot
    if all_pass
        println("  ALL CHECKS PASSED within $(tol) tolerance.\n")
    else
        println("  SOME CHECKS FAILED. Investigate numerical differences.\n")

        # Print first few differing points for debugging
        ref_I_tol_arr = ref_I_tol
        for (name, impl, ref_vals) in [("I_tol", I_tol, ref_I_tol),
                                         ("I_seq", I_seq, ref_I_seq),
                                         ("I_cot", I_cot, ref_I_cot)]
            diff = abs.(impl .- ref_vals)
            denom = max.(abs.(ref_vals), 1e-30)
            rel_err = diff ./ denom
            bad_idx = findall(rel_err .>= tol)
            if !isempty(bad_idx)
                println("  $name: $(length(bad_idx)) points exceed tolerance")
                for idx in bad_idx[1:min(5, length(bad_idx))]
                    @printf("    Vsd[%d]=%.4f: julia=%.15e  matlab=%.15e  relerr=%.4e\n",
                            idx, Vsd[idx], impl[idx], ref_vals[idx], rel_err[idx])
                end
            end
        end
    end
end

# Run
main()
