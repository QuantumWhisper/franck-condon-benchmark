using SpecialFunctions: digamma, trigamma

# ============================================================================
# High-precision digamma/trigamma via BigFloat asymptotic series.
# This matches MATLAB's approach of using symbolic math (arbitrary precision)
# for the polygamma function. The regularized integrals involve cancellation
# of digamma values, amplifying any precision loss at Float64. Using BigFloat
# internally eliminates this amplification.
# ============================================================================

# First 20 even Bernoulli numbers B_2, B_4, ..., B_40 (precomputed as BigFloat)
const _BERNOULLI_EVEN = let
    # Exact rational values of B_{2k} for k=1..20
    rats = [1//6, -1//30, 1//42, -1//30, 5//66, -691//2730, 7//6, -3617//510,
            43867//798, -174611//330, 854513//138, -236364091//2730, 8553103//6,
            -23749461029//870, 8615841276005//14322, -7709321041217//510,
            2577687858367//6, -26315271553053477373//1919190,
            2929993913841559//6, -261082718496449122051//13530]
    BigFloat.(rats)
end

"""
    _digamma_bf(z::Complex{BigFloat}) -> Complex{BigFloat}

Digamma function for complex BigFloat via asymptotic series with recurrence shift.
ψ(z) = ln(z) - 1/(2z) - Σ_{k=1}^{K} B_{2k}/(2k·z^{2k})
"""
function _digamma_bf(z::Complex{BigFloat})::Complex{BigFloat}
    # Reflection formula for Re(z) <= 0
    if real(z) <= 0
        # ψ(z) = ψ(1-z) - π/tan(πz)
        return _digamma_bf(one(z) - z) - big(π) * cos(big(π) * z) / sin(big(π) * z)
    end

    # Recurrence shift: ψ(z) = ψ(z+1) - 1/z until |z| >= 20
    result = zero(z)
    z_s = z
    while abs(z_s) < 20
        result -= 1 / z_s
        z_s += 1
    end

    # Asymptotic expansion
    result += log(z_s) - 1 / (2 * z_s)
    z_sq = z_s * z_s
    power = z_sq
    for (k, bk) in enumerate(_BERNOULLI_EVEN)
        result -= bk / (2k * power)
        power *= z_sq
    end

    return result
end

"""
    _trigamma_bf(z::Complex{BigFloat}) -> Complex{BigFloat}

Trigamma function for complex BigFloat via asymptotic series.
ψ'(z) = 1/z + 1/(2z²) + Σ_{k=1}^{K} B_{2k}/z^{2k+1}
"""
function _trigamma_bf(z::Complex{BigFloat})::Complex{BigFloat}
    # Reflection formula for Re(z) <= 0
    if real(z) <= 0
        # ψ'(1-z) = (π/sin(πz))² - ψ'(z)  →  ψ'(z) = (π/sin(πz))² - ψ'(1-z)
        return (big(π) / sin(big(π) * z))^2 - _trigamma_bf(one(z) - z)
    end

    # Recurrence shift: ψ'(z) = ψ'(z+1) + 1/z² until |z| >= 20
    result = zero(z)
    z_s = z
    while abs(z_s) < 20
        result += 1 / (z_s * z_s)
        z_s += 1
    end

    # Asymptotic expansion
    result += 1 / z_s + 1 / (2 * z_s * z_s)
    z_sq = z_s * z_s
    power = z_sq * z_s  # z^3
    for bk in _BERNOULLI_EVEN
        result += bk / power
        power *= z_sq
    end

    return result
end

"""
    digamma_hp(z::Complex{Float64}) -> Complex{Float64}

High-precision digamma: compute in BigFloat, return Float64.
Matches MATLAB's `double(psi(0, sym(z)))`.
"""
function digamma_hp(z::Complex{Float64})::Complex{Float64}
    Complex{Float64}(_digamma_bf(Complex{BigFloat}(z)))
end

"""
    trigamma_hp(z::Complex{Float64}) -> Complex{Float64}

High-precision trigamma: compute in BigFloat, return Float64.
Matches MATLAB's `double(psi(1, sym(z)))`.
"""
function trigamma_hp(z::Complex{Float64})::Complex{Float64}
    Complex{Float64}(_trigamma_bf(Complex{BigFloat}(z)))
end

"""
    regularized_I(E1, E2, epsilon1_in, epsilon2_in, T) -> Matrix{Float64}

Analytically regularized cotunneling integral I.
Matches MATLAB `regularizedI.m` exactly, including broadcasting behavior.

epsilon1 is transposed internally (becomes column), epsilon2 stays as row.
Returns a matrix of size [length(epsilon1), length(epsilon2)].
"""
function regularized_I(E1, E2, epsilon1_in::AbstractVector, epsilon2_in::AbstractVector, T)
    beta = 1.0 / (KB_EV * T)

    # MATLAB: epsilon1 = epsilon1' (transpose to column for broadcasting)
    eps1 = reshape(collect(Float64, epsilon1_in), :, 1)  # column [N, 1]
    eps2 = reshape(collect(Float64, epsilon2_in), 1, :)   # row [1, N]

    a1 = @. 0.5 + im * beta * (E2 - eps1) / (2π)
    a2 = @. 0.5 - im * beta * (E2 - eps2) / (2π)
    a3 = @. 0.5 + im * beta * (E1 - eps1) / (2π)
    a4 = @. 0.5 - im * beta * (E1 - eps2) / (2π)

    t1 = digamma.(a1)
    t2 = digamma.(a2)
    t3 = digamma.(a3)
    t4 = digamma.(a4)

    bose_val = bose_fcn(E2 - E1, T)

    return @. bose_val / (eps1 - eps2) * real(t1 - t2 - t3 + t4)
end

"""
    regularized_J(E1, E2, epsilon_in, T) -> Matrix or Vector

Analytically regularized cotunneling integral J.
Matches MATLAB `regularizedJ.m` exactly, including broadcasting behavior.

epsilon is transposed internally (column or transposed matrix).
E2 can be scalar or vector (reshaped to row for broadcasting).
"""
function regularized_J(E1, E2_in, epsilon_in, T)
    beta = 1.0 / (KB_EV * T)

    # MATLAB: epsilon = epsilon' (transpose)
    if epsilon_in isa AbstractMatrix
        epsilon = permutedims(epsilon_in)  # transpose matrix
    else
        epsilon = reshape(collect(Float64, epsilon_in), :, 1)  # vector → column
    end

    # E2: ensure row shape for broadcasting
    if E2_in isa AbstractVector
        E2 = reshape(collect(Float64, E2_in), 1, :)
    else
        E2 = E2_in  # scalar
    end

    a1 = @. 0.5 + im * beta * (E2 - epsilon) / (2π)
    a2 = @. 0.5 + im * beta * (E1 - epsilon) / (2π)

    t1 = trigamma.(a1)
    t2 = trigamma.(a2)

    bose_val = bose_fcn.(E2 .- E1, T)

    return @. beta / (2π) * bose_val * imag(t1 - t2)
end
