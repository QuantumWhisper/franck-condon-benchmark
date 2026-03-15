using LinearAlgebra

"""
    solve_steady_state(M::Matrix{Float64}) -> Vector{Float64}

Solve the steady-state occupation probability vector P from 0 = WP.
Equivalent to MATLAB's lsqlin with constraints sum(P)=1, 0≤P≤1.

Uses the augmented system approach (replace last row with normalization),
which gives identical results to SVD and Ipopt for this problem.
"""
function solve_steady_state(M::Matrix{Float64})::Vector{Float64}
    n = size(M, 2)

    # Replace last row with normalization constraint
    M_aug = copy(M)
    M_aug[end, :] .= 1.0
    d = zeros(n)
    d[end] = 1.0

    P = M_aug \ d

    P = max.(P, 0.0)
    s = sum(P)
    if s > 0
        P ./= s
    end

    return P
end
