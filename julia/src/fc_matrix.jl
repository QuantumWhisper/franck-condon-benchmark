using SpecialFunctions: gamma

"""
    fc_matrix_single(q1::Int, q2::Int, lambda::Float64) -> Float64

Single Franck-Condon matrix element <q2|D(lambda)|q1>.
Matches MATLAB `FCMatrixSingle.m` exactly.

Formula: sign(q2-q1)^(q1-q2) * lambda^(Q-q) * exp(-lambda^2/2) * sqrt(q!/Q!) * L_q^{Q-q}(lambda^2)
where q = min(q1,q2), Q = max(q1,q2).
"""
function fc_matrix_single(q1::Int, q2::Int, lambda::Float64)::Float64
    q = min(q1, q2)
    Q = max(q1, q2)

    l = laguerre_L(q, Q - q, lambda^2)

    # sign(q2-q1)^(q1-q2):
    #   q1==q2: 0^0 = 1
    #   q2>q1:  1^(negative) = 1
    #   q2<q1:  (-1)^(positive integer)
    if q1 == q2
        sf = 1.0
    elseif q2 > q1
        sf = 1.0
    else  # q2 < q1
        sf = iseven(q1 - q2) ? 1.0 : -1.0
    end

    # Use gamma(n+1) = n! to match MATLAB factorial behavior (Inf for large n)
    fq = gamma(Float64(q + 1))
    fQ = gamma(Float64(Q + 1))

    M = sf * lambda^(Q - q) * exp(-lambda^2 / 2) * sqrt(fq / fQ) * l

    # Handle NaN/Inf from overflow (MATLAB returns 0 in these cases)
    return (isnan(M) || isinf(M)) ? 0.0 : M
end

const FC_MAX_N = 256

"""
FC matrix cache for a fixed lambda value. Pre-populated flat Matrix for O(1) lookup.
Thread-safe after construction: the cache Matrix is populated once and then only read.
"""
struct FCCache
    cache::Matrix{Float64}   # FC_MAX_N × FC_MAX_N, indexed as [q1+1, q2+1]
    lambda::Float64
end

"""
    FCCache(lambda) -> FCCache

Construct and pre-populate the FC cache for all q1, q2 in 0:FC_MAX_N-1.
Cost: ~65536 evaluations of fc_matrix_single (~5-50 ms). One-time cost per simulation.

FC_MAX_N=256 covers ALL non-zero FC matrix elements: gamma(172) overflows in Float64,
so fc_matrix_single returns 0.0 for any (q1,q2) where max(q1,q2) ≥ 171. The convergence
wrappers may access q ≥ 256, but those entries are guaranteed 0.0 (handled by get_fc!).

After construction, the cache is immutable and thread-safe for concurrent reads.
"""
function FCCache(lambda::Float64)
    cache = Matrix{Float64}(undef, FC_MAX_N, FC_MAX_N)
    for q2 in 0:FC_MAX_N-1
        for q1 in 0:FC_MAX_N-1
            cache[q1 + 1, q2 + 1] = fc_matrix_single(q1, q2, lambda)
        end
    end
    return FCCache(cache, lambda)
end

"""
    get_fc!(fc::FCCache, q1::Int, q2::Int) -> Float64

Retrieve pre-computed FC matrix element. O(1) array access, no allocation.
Name kept with `!` for backward compatibility; function is actually read-only
after FCCache construction.

For q ≥ FC_MAX_N (256), returns 0.0 directly. This is mathematically exact:
gamma(172) overflows in Float64, so fc_matrix_single returns 0.0 for any
(q1, q2) where max(q1, q2) ≥ 171. Since FC_MAX_N=256 > 171, all non-zero
FC elements are in the cache.
"""
@inline function get_fc!(fc::FCCache, q1::Int, q2::Int)::Float64
    if q1 >= FC_MAX_N || q2 >= FC_MAX_N
        return 0.0
    end
    @inbounds return fc.cache[q1 + 1, q2 + 1]
end

"""
    fc_matrix(q1_range, q2_range, fc::FCCache) -> Matrix{Float64}

Compute FC matrix M[i,j] = FC(q1_range[i], q2_range[j], lambda).
Matches MATLAB `FCMatrix(q1, q2, lambda)`.
"""
function fc_matrix(q1_range, q2_range, fc::FCCache)::Matrix{Float64}
    nq1 = length(q1_range)
    nq2 = length(q2_range)
    M = Matrix{Float64}(undef, nq1, nq2)
    @inbounds for (j, q2) in enumerate(q2_range)
        for (i, q1) in enumerate(q1_range)
            M[i, j] = get_fc!(fc, q1, q2)
        end
    end
    return M
end
