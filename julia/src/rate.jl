"""
Rate store: pre-computed rates indexed by (n1, n2, q1, lead) → Vector{Float64} of length N.
"""
const RateStore = Dict{Tuple{Int,Int,Int,Int},Vector{Float64}}

"""
    spin_degeneracy(n1, n2) -> Float64

Spin degeneracy factor s(n1, n2). Matches MATLAB `s_` in m_rateW.m.
"""
function spin_degeneracy(n1::Int, n2::Int)::Float64
    if n1 == 0
        return 2.0
    else  # n1 == 1
        return n2 == 0 ? 1.0 : 2.0
    end
end

"""
    m_rateW(n1, n2, q1, q2_vec, vmode, alphaL, alphaR, lambda, Vsd, T, eta, lead, Vg, fc)

Core rate computation. Matches MATLAB `m_rateW.m` exactly.
Returns vector of rates for each q2 value.
"""
function m_rateW(n1::Int, n2::Int, q1::Int, q2_vec::Vector{Int},
                 vmode::Float64, alphaL::Float64, alphaR::Float64, lambda::Float64,
                 Vsd::Float64, T::Float64, eta::Float64, lead::Int, Vg::Float64,
                 fc::FCCache)
    gammaL = alphaL * vmode
    gammaR = alphaR * vmode
    epsilond = 0.0 + Vg
    muL = eta * Vsd
    muR = -(1 - eta) * Vsd

    s = spin_degeneracy(n1, n2)
    gamma = lead == 1 ? gammaL : gammaR
    mu = lead == 1 ? muL : muR

    local w::Vector{Float64}

    if n1 == 1 && n2 == 0
        # Sequential tunneling: 1→0 (occupied to empty)
        fc_vals = Float64[get_fc!(fc, q1, q) for q in q2_vec]
        w = s * gamma / HBAR_EV .* fc_vals.^2 .*
            (1.0 .- fermi.(epsilond .- (q2_vec .- q1) .* vmode, mu, T))

    elseif n1 == 0 && n2 == 1
        # Sequential tunneling: 0→1 (empty to occupied)
        fc_vals = Float64[get_fc!(fc, q1, q) for q in q2_vec]
        w = s * gamma / HBAR_EV .* fc_vals.^2 .*
            fermi.(epsilond .+ (q2_vec .- q1) .* vmode, mu, T)

    elseif n1 == 0 && n2 == 0
        # Cotunneling: 0→0
        if lead == 1  # L → R
            w = s / (2π * HBAR_EV) * gammaL * gammaR .*
                (sumMMr(q1, q2_vec, lambda, muL, muR, vmode, epsilond, T, fc) .+
                 sumMMMMrs(q1, q2_vec, lambda, muL, muR, vmode, epsilond, T, fc))
        else  # R → L
            w = s / (2π * HBAR_EV) * gammaL * gammaR .*
                (sumMMr(q1, q2_vec, lambda, muR, muL, vmode, epsilond, T, fc) .+
                 sumMMMMrs(q1, q2_vec, lambda, muR, muL, vmode, epsilond, T, fc))
        end

    else  # n1 == 1 && n2 == 1
        # Cotunneling: 1→1
        if lead == 1
            w = s / (2π * HBAR_EV) * gammaL * gammaR .*
                (sumMMr11(q1, q2_vec, lambda, muL, muR, vmode, epsilond, T, fc) .+
                 sumMMMMrs11(q1, q2_vec, lambda, muL, muR, vmode, epsilond, T, fc))
        else
            w = s / (2π * HBAR_EV) * gammaL * gammaR .*
                (sumMMr11(q1, q2_vec, lambda, muR, muL, vmode, epsilond, T, fc) .+
                 sumMMMMrs11(q1, q2_vec, lambda, muR, muL, vmode, epsilond, T, fc))
        end
    end

    # Replace NaN with 0 (matches MATLAB: w(isnan(w)) = 0)
    replace!(x -> isnan(x) ? 0.0 : x, w)
    return w
end

"""
    calculate_all_rateW!(store, N, vmode, alphaL, alphaR, lambda, Vsd, T, eta, lead, Vg, fc)

Pre-compute all rates for a given lead and bias point. Fills the rate store.
Matches MATLAB `calculateAllRateW.m`.
"""
function calculate_all_rateW!(store::RateStore, N::Int, vmode::Float64, alphaL::Float64,
                               alphaR::Float64, lambda::Float64, Vsd::Float64, T::Float64,
                               eta::Float64, lead::Int, Vg::Float64, fc::FCCache)
    q2_vec = collect(0:N-1)
    for n1 in 0:1
        for n2 in 0:1
            for ii in 1:N
                q1 = ii - 1
                w = m_rateW(n1, n2, q1, q2_vec, vmode, alphaL, alphaR, lambda,
                            Vsd, T, eta, lead, Vg, fc)
                store[(n1, n2, q1, lead)] = w
            end
        end
    end
end

"""
    rateW_from_store(store, n1, n2, q1, q2_indices, lead)

Retrieve pre-computed rates from store. q2_indices are 0-based.
"""
function rateW_from_store(store::RateStore, n1::Int, n2::Int, q1::Int,
                          q2_indices, lead::Int)
    w_full = store[(n1, n2, q1, lead)]
    if q2_indices isa Integer
        return w_full[q2_indices + 1]
    else
        return w_full[collect(q2_indices) .+ 1]
    end
end

"""
    rateW_lead_from_store(store, n1, n2, q1, q2)

Rate summed over both leads. q2 can be scalar or vector (0-based).
Matches MATLAB `rateW_lead.m`.
"""
function rateW_lead_from_store(store::RateStore, n1::Int, n2::Int, q1::Int, q2)
    wl = rateW_from_store(store, n1, n2, q1, q2, 1)
    wr = rateW_from_store(store, n1, n2, q1, q2, -1)
    return wl .+ wr
end
