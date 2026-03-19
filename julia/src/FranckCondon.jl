module FranckCondon

using SpecialFunctions
using LinearAlgebra
using Printf

include("constants.jl")
include("laguerre.jl")
include("fc_matrix.jl")
include("fermi_bose.jl")
include("regularized.jl")
include("cotunneling.jl")
include("rate.jl")
include("matrix.jl")
include("solver.jl")
include("current.jl")

export simulate_iv, FCCache, ELEMENTARY_CHARGE

"""
    simulate_iv(N, vmode, alphaL, alphaR, lambda, Vsd_vec, T, eta, Vg, tau; verbose=true)

Run the full I-V simulation. Returns (Vsd_vec, I_tol, I_seq, I_cot) in SI units (Amperes).
Matches MATLAB `run_benchmark.m` simulation loop exactly.

Threading: when Julia is started with multiple threads (julia --threads=auto),
the bias-point loop runs in parallel. The FCCache is pre-populated once and shared
read-only across all threads. Each thread gets its own RateStore (no synchronization).
"""
function simulate_iv(N::Int, vmode::Float64, alphaL::Float64, alphaR::Float64,
                     lambda::Float64, Vsd_vec::Vector{Float64}, T::Float64, eta::Float64,
                     Vg::Float64, tau::Float64; verbose::Bool=true)
    nVsd = length(Vsd_vec)
    I_tol = zeros(nVsd)
    I_seq = zeros(nVsd)
    I_cot = zeros(nVsd)

    # FC matrix cache: pre-populated in constructor, shared read-only across threads
    fc = FCCache(lambda)

    nthreads = Threads.nthreads()

    if verbose && nthreads > 1
        @printf("Parallel: %d threads × %d bias points\n", nthreads, nVsd)
    end

    Threads.@threads for vv in 1:nVsd
        v = Vsd_vec[vv]

        # Per-bias-point verbose output (only when single-threaded to avoid interleaving)
        if verbose && nthreads == 1
            @printf("Vsd = %.4f V  (%d/%d)\n", v, vv, nVsd)
        end

        # Per-iteration rate store (each bias point is independent)
        store = RateStore(N)
        calculate_all_rateW!(store, N, vmode, alphaL, alphaR, lambda, v, T, eta, 1, Vg, fc)
        calculate_all_rateW!(store, N, vmode, alphaL, alphaR, lambda, v, T, eta, -1, Vg, fc)

        # Compute current
        i_tol, i_seq, i_cot = current_from_rate_equations(
            v, N, vmode, alphaL, alphaR, lambda, T, eta, Vg, tau, store)

        # Sign convention: current flows opposite to bias
        if v != 0.0
            s = -sign(v)
            i_tol *= s
            i_seq *= s
            i_cot *= s
        end

        # Write to output arrays (each vv is unique — no data race)
        I_tol[vv] = i_tol
        I_seq[vv] = i_seq
        I_cot[vv] = i_cot
    end

    # Convert to SI (Amperes)
    e = ELEMENTARY_CHARGE
    I_tol .*= e
    I_seq .*= e
    I_cot .*= e

    return (Vsd_vec, I_tol, I_seq, I_cot)
end

end # module
