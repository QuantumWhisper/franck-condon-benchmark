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
"""
function simulate_iv(N::Int, vmode::Float64, alphaL::Float64, alphaR::Float64,
                     lambda::Float64, Vsd_vec::Vector{Float64}, T::Float64, eta::Float64,
                     Vg::Float64, tau::Float64; verbose::Bool=true)
    nVsd = length(Vsd_vec)
    I_tol = zeros(nVsd)
    I_seq = zeros(nVsd)
    I_cot = zeros(nVsd)

    # FC matrix cache persists across all bias points
    fc = FCCache(lambda)

    for vv in 1:nVsd
        v = Vsd_vec[vv]
        if verbose
            @printf("Vsd = %.4f V  (%d/%d)\n", v, vv, nVsd)
        end

        # Pre-compute rate matrices for both leads
        store = RateStore()
        calculate_all_rateW!(store, N, vmode, alphaL, alphaR, lambda, v, T, eta, 1, Vg, fc)
        calculate_all_rateW!(store, N, vmode, alphaL, alphaR, lambda, v, T, eta, -1, Vg, fc)

        # Compute current
        I_tol[vv], I_seq[vv], I_cot[vv] = current_from_rate_equations(
            v, N, vmode, alphaL, alphaR, lambda, T, eta, Vg, tau, store)

        # Sign convention: current flows opposite to bias
        if v != 0.0
            s = -sign(v)
            I_tol[vv] *= s
            I_seq[vv] *= s
            I_cot[vv] *= s
        end
    end

    # Convert to SI (Amperes)
    e = ELEMENTARY_CHARGE
    I_tol .*= e
    I_seq .*= e
    I_cot .*= e

    return (Vsd_vec, I_tol, I_seq, I_cot)
end

end # module
