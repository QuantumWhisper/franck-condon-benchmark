"""
    laguerre_L(n::Int, alpha::Int, x::Float64) -> Float64

Generalized (associated) Laguerre polynomial L_n^alpha(x) via the standard three-term
forward recurrence. Matches MATLAB's `laguerreL(n, alpha, x)` to full Float64 precision.

    L_0^alpha(x) = 1
    L_1^alpha(x) = 1 + alpha - x
    L_{k+1}^alpha(x) = ((2k + 1 + alpha - x) * L_k - (k + alpha) * L_{k-1}) / (k + 1)
"""
function laguerre_L(n::Int, alpha::Int, x::Float64)::Float64
    n == 0 && return 1.0
    L_prev = 1.0
    L_curr = 1.0 + alpha - x
    n == 1 && return L_curr
    for k in 1:(n - 1)
        L_next = ((2k + 1 + alpha - x) * L_curr - (k + alpha) * L_prev) / (k + 1)
        L_prev = L_curr
        L_curr = L_next
    end
    return L_curr
end
