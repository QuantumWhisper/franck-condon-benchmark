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

"""
FC matrix cache for a fixed lambda value. Avoids recomputation of FC elements.
"""
mutable struct FCCache
    cache::Dict{Tuple{Int,Int},Float64}
    lambda::Float64
end

FCCache(lambda::Float64) = FCCache(Dict{Tuple{Int,Int},Float64}(), lambda)

function get_fc!(fc::FCCache, q1::Int, q2::Int)::Float64
    get!(fc.cache, (q1, q2)) do
        fc_matrix_single(q1, q2, fc.lambda)
    end
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
    for (j, q2) in enumerate(q2_range)
        for (i, q1) in enumerate(q1_range)
            M[i, j] = get_fc!(fc, q1, q2)
        end
    end
    return M
end
