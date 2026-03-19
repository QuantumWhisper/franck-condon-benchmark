"""
Rate store: pre-computed rates in a flat array indexed by (n1, n2, q1, lead, q2).
Replaces Dict{Tuple{Int,Int,Int,Int},Vector{Float64}} for O(1) cache-friendly access.

Index formula (matching Rust's rate.rs):
  flat_idx = (((n1*2 + n2)*N + q1)*2 + lead_idx)*N + q2 + 1
  where lead_idx = 0 for lead=+1, 1 for lead=-1.
Total size: 2 × 2 × N × 2 × N = 8N².
"""
struct RateStore
    n::Int
    data::Vector{Float64}
end

"""
    RateStore(n::Int) -> RateStore

Create a rate store for N phonon states. Allocates flat array of size 8N².
"""
RateStore(n::Int) = RateStore(n, zeros(Float64, 2 * 2 * n * 2 * n))

@inline _lead_idx(lead::Int) = lead == 1 ? 0 : 1

@inline function _rate_idx(n::Int, n1::Int, n2::Int, q1::Int, lead_idx::Int, q2::Int)::Int
    (((n1 * 2 + n2) * n + q1) * 2 + lead_idx) * n + q2 + 1
end

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
    lid = _lead_idx(lead)
    n = store.n
    for n1 in 0:1
        for n2 in 0:1
            for ii in 1:N
                q1 = ii - 1
                w = m_rateW(n1, n2, q1, q2_vec, vmode, alphaL, alphaR, lambda,
                            Vsd, T, eta, lead, Vg, fc)
                # Store rates in flat array (contiguous block for all q2)
                base = _rate_idx(n, n1, n2, q1, lid, 0)
                @inbounds for i in 1:n
                    store.data[base + i - 1] = w[i]
                end
            end
        end
    end
end

"""
    rateW_from_store(store, n1, n2, q1, q2, lead) -> Float64

Retrieve single pre-computed rate from store. q2 is 0-based.
"""
@inline function rateW_from_store(store::RateStore, n1::Int, n2::Int, q1::Int,
                                   q2::Int, lead::Int)::Float64
    @inbounds return store.data[_rate_idx(store.n, n1, n2, q1, _lead_idx(lead), q2)]
end

"""
    rateW_from_store(store, n1, n2, q1, q2_vec, lead) -> SubArray

Retrieve pre-computed rates for all q2 values. Returns a view (zero allocation).
q2_vec must be 0:N-1 (the standard call pattern).
"""
@inline function rateW_from_store(store::RateStore, n1::Int, n2::Int, q1::Int,
                                   q2_vec::AbstractVector{Int}, lead::Int)
    n = store.n
    base = _rate_idx(n, n1, n2, q1, _lead_idx(lead), 0)
    @inbounds return @view store.data[base:base+n-1]
end

"""
    rateW_lead_from_store(store, n1, n2, q1, q2)

Rate summed over both leads. q2 can be scalar (0-based) or vector.
Matches MATLAB `rateW_lead.m`.
"""
@inline function rateW_lead_from_store(store::RateStore, n1::Int, n2::Int, q1::Int, q2::Int)
    wl = rateW_from_store(store, n1, n2, q1, q2, 1)
    wr = rateW_from_store(store, n1, n2, q1, q2, -1)
    return wl + wr
end

function rateW_lead_from_store(store::RateStore, n1::Int, n2::Int, q1::Int, q2_vec::AbstractVector{Int})
    wl = rateW_from_store(store, n1, n2, q1, q2_vec, 1)
    wr = rateW_from_store(store, n1, n2, q1, q2_vec, -1)
    return wl .+ wr
end
