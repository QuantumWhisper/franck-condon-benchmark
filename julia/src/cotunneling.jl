# ==================== Inner computation functions ====================

"""
    m_sumMMr(N, q1, q2_vec, lambda, muL, muR, vmode, epsilond, T, fc)

Single cotunneling sum over virtual states r for n=0→0 transition.
Matches MATLAB `m_sumMMr.m` exactly.
Returns vector of length nq2.
"""
function m_sumMMr(N::Int, q1::Int, q2_vec::Vector{Int}, lambda::Float64,
                  muL::Float64, muR::Float64, vmode::Float64, epsilond::Float64,
                  T::Float64, fc::FCCache)
    rr = collect(0:N-1)
    nq2 = length(q2_vec)

    # FCMatrix(q2, rr)' .* FCMatrix(q1, rr)' → [nr, nq2]
    fc_q2_rr = fc_matrix(q2_vec, rr, fc)     # [nq2, nr]
    fc_q1_rr = fc_matrix([q1], rr, fc)        # [1, nr]
    # MATLAB: abs(FCMatrix(q2,rr,lambda)'.*FCMatrix(q1,rr,lambda)').^2
    MM_square = abs.(permutedims(fc_q2_rr) .* permutedims(fc_q1_rr)).^2  # [nr, nq2]

    # Jr = regularizedJ(muL, muR-(q1-q2)*vmode, epsilond-(q1-rr)*vmode, T)
    E2_vec = muR .- (q1 .- q2_vec) .* vmode   # [nq2]
    eps_vec = epsilond .- (q1 .- rr) .* vmode  # [nr]
    Jr = regularized_J(muL, E2_vec, eps_vec, T)  # [nr, nq2]

    target = MM_square .* Jr
    return vec(sum(target, dims=1))  # [nq2]
end

"""
    m_sumMMr11(N, q1, q2_vec, lambda, muL, muR, vmode, epsilond, T, fc)

Single cotunneling sum for n=1→1 transition.
Matches MATLAB `m_sumMMr11.m`.
"""
function m_sumMMr11(N::Int, q1::Int, q2_vec::Vector{Int}, lambda::Float64,
                    muL::Float64, muR::Float64, vmode::Float64, epsilond::Float64,
                    T::Float64, fc::FCCache)
    rr = collect(0:N-1)
    nq2 = length(q2_vec)

    # MATLAB: abs(FCMatrix(q2,rr).*FCMatrix(q1,rr)).^2'
    fc_q2_rr = fc_matrix(q2_vec, rr, fc)   # [nq2, nr]
    fc_q1_rr = fc_matrix([q1], rr, fc)      # [1, nr]
    MM_square = permutedims(abs.(fc_q2_rr .* fc_q1_rr).^2)  # [nr, nq2]

    # epsilon = epsilond + (q2' - rr) * vmode → matrix [nq2, nr]
    q2_col = reshape(q2_vec, :, 1)
    rr_row = reshape(rr, 1, :)
    eps_matrix = epsilond .+ (q2_col .- rr_row) .* vmode  # [nq2, nr]

    E2_vec = muR .- (q1 .- q2_vec) .* vmode
    Jr = regularized_J(muL, E2_vec, eps_matrix, T)  # [nr, nq2] after internal transpose

    target = MM_square .* Jr
    return vec(sum(target, dims=1))
end

"""
    m_sumMMMMrs(N, q1, q2, lambda, muL, muR, vmode, epsilond, T, fc)

Double cotunneling sum over r,s for n=0→0 transition.
Matches MATLAB `m_sumMMMMrs.m` + `m_MMMM_quad_Irs.m`.
Returns scalar.
"""
function m_sumMMMMrs(N::Int, q1::Int, q2::Int, lambda::Float64,
                     muL::Float64, muR::Float64, vmode::Float64, epsilond::Float64,
                     T::Float64, fc::FCCache)
    rr = collect(0:N-1)
    ss = collect(0:N-1)
    nr = N

    # MMMM_quad2D = FC(q2,rr)' .* conj(FC(q1,rr))' .* conj(FC(q2,ss)) .* FC(q1,ss)
    fc_q2_rr = fc_matrix([q2], rr, fc)  # [1, nr]
    fc_q1_rr = fc_matrix([q1], rr, fc)  # [1, nr]
    fc_q2_ss = fc_matrix([q2], ss, fc)  # [1, nr]
    fc_q1_ss = fc_matrix([q1], ss, fc)  # [1, nr]

    MMMM = permutedims(fc_q2_rr) .* permutedims(conj.(fc_q1_rr)) .*
           conj.(fc_q2_ss) .* fc_q1_ss  # [nr, nr]

    # Irs = regularizedI(muL, muR-(q1-q2)*vmode, epsilond-(q1-rr)*vmode, epsilond-(q1-ss)*vmode, T)
    E2 = muR - (q1 - q2) * vmode
    eps1 = epsilond .- (q1 .- rr) .* vmode
    eps2 = epsilond .- (q1 .- ss) .* vmode
    Irs = regularized_I(muL, E2, eps1, eps2, T)  # [nr, nr]

    target = MMMM .* Irs
    # Zero diagonal
    @inbounds for i in 1:nr
        target[i, i] = 0.0
    end
    return sum(target)
end

"""
    m_sumMMMMrs11(N, q1, q2, lambda, muL, muR, vmode, epsilond, T, fc)

Double cotunneling sum for n=1→1 transition.
Matches MATLAB `m_sumMMMMrs11.m`.
"""
function m_sumMMMMrs11(N::Int, q1::Int, q2::Int, lambda::Float64,
                       muL::Float64, muR::Float64, vmode::Float64, epsilond::Float64,
                       T::Float64, fc::FCCache)
    rr = collect(0:N-1)
    ss = collect(0:N-1)
    nr = N

    fc_q2_rr = fc_matrix([q2], rr, fc)
    fc_q1_rr = fc_matrix([q1], rr, fc)
    fc_q2_ss = fc_matrix([q2], ss, fc)
    fc_q1_ss = fc_matrix([q1], ss, fc)

    MMMM = permutedims(fc_q2_rr) .* permutedims(conj.(fc_q1_rr)) .*
           conj.(fc_q2_ss) .* fc_q1_ss

    # Different epsilon for n=1→1: epsilond+(q2-rr)*vmode instead of epsilond-(q1-rr)*vmode
    E2 = muR - (q1 - q2) * vmode
    eps1 = epsilond .+ (q2 .- rr) .* vmode
    eps2 = epsilond .+ (q2 .- ss) .* vmode
    Irs = regularized_I(muL, E2, eps1, eps2, T)

    target = MMMM .* Irs
    @inbounds for i in 1:nr
        target[i, i] = 0.0
    end
    return sum(target)
end

# ==================== Convergence wrappers ====================

"""
    sumMMr(q1, q2_vec, lambda, muL, muR, vmode, epsilond, T, fc)

Converged single cotunneling sum for n=0→0. Uses adaptive N with log10 convergence criterion.
Matches MATLAB `sumMMr.m`.
"""
function sumMMr(q1::Int, q2_vec::Vector{Int}, lambda::Float64,
                muL::Float64, muR::Float64, vmode::Float64, epsilond::Float64,
                T::Float64, fc::FCCache)
    tempN = round(Int, lambda^2.2 * 3)
    epsilon = 1e-14

    temp_tol = m_sumMMr(tempN, q1, q2_vec, lambda, muL, muR, vmode, epsilond, T, fc)
    tempN += 5
    temp_tol2 = m_sumMMr(tempN, q1, q2_vec, lambda, muL, muR, vmode, epsilond, T, fc)

    # MATLAB uses log10 for convergence (not log)
    relative_diff = abs.(log10.(complex.(temp_tol2 ./ temp_tol)))

    while any(x -> isfinite(x) && x > epsilon, relative_diff)
        locs = [isfinite(rd) && rd > epsilon for rd in relative_diff]
        tempN += min(max(round(Int, tempN * 0.5), 10), 20)
        temp_tol .= temp_tol2
        q2_subset = q2_vec[locs]
        temp_tol2_partial = m_sumMMr(tempN, q1, q2_subset, lambda, muL, muR, vmode, epsilond, T, fc)
        temp_tol2[locs] .= temp_tol2_partial
        relative_diff = abs.(log10.(complex.(temp_tol2 ./ temp_tol)))
    end

    return temp_tol2
end

"""
    sumMMr11(q1, q2_vec, lambda, muL, muR, vmode, epsilond, T, fc)

Converged single cotunneling sum for n=1→1. Matches MATLAB `sumMMr11.m`.
"""
function sumMMr11(q1::Int, q2_vec::Vector{Int}, lambda::Float64,
                  muL::Float64, muR::Float64, vmode::Float64, epsilond::Float64,
                  T::Float64, fc::FCCache)
    tempN = round(Int, lambda^2.2 * 3)
    epsilon = 1e-14

    temp_tol = m_sumMMr11(tempN, q1, q2_vec, lambda, muL, muR, vmode, epsilond, T, fc)
    tempN += 5
    temp_tol2 = m_sumMMr11(tempN, q1, q2_vec, lambda, muL, muR, vmode, epsilond, T, fc)

    relative_diff = abs.(log10.(complex.(temp_tol2 ./ temp_tol)))

    while any(x -> isfinite(x) && x > epsilon, relative_diff)
        locs = [isfinite(rd) && rd > epsilon for rd in relative_diff]
        tempN += min(max(round(Int, tempN * 0.5), 10), 20)
        temp_tol .= temp_tol2
        q2_subset = q2_vec[locs]
        temp_tol2_partial = m_sumMMr11(tempN, q1, q2_subset, lambda, muL, muR, vmode, epsilond, T, fc)
        temp_tol2[locs] .= temp_tol2_partial
        relative_diff = abs.(log10.(complex.(temp_tol2 ./ temp_tol)))
    end

    return temp_tol2
end

"""
    sumMMMMrs(q1, q2_vec, lambda, muL, muR, vmode, epsilond, T, fc)

Converged double cotunneling sum for n=0→0.
Processes each q2 individually. Uses log (natural) convergence criterion.
Matches MATLAB `sumMMMMrs.m`.
"""
function sumMMMMrs(q1::Int, q2_vec::Vector{Int}, lambda::Float64,
                   muL::Float64, muR::Float64, vmode::Float64, epsilond::Float64,
                   T::Float64, fc::FCCache)
    nq2 = length(q2_vec)
    result = zeros(nq2)

    for qq2 in 1:nq2
        q2 = q2_vec[qq2]
        tempN = round(Int, lambda^2 * 4)
        eps_conv = 1e-14

        temp_tol = m_sumMMMMrs(tempN, q1, q2, lambda, muL, muR, vmode, epsilond, T, fc)
        tempN += 5
        temp_tol2 = m_sumMMMMrs(tempN, q1, q2, lambda, muL, muR, vmode, epsilond, T, fc)

        # MATLAB uses log (natural) for convergence
        rd = temp_tol == 0.0 && temp_tol2 == 0.0 ? 0.0 :
             temp_tol == 0.0 ? Inf :
             abs(log(complex(temp_tol2 / temp_tol)))

        while isfinite(rd) && rd > eps_conv
            temp_tol = temp_tol2
            tempN += min(max(round(Int, tempN * 0.5), 20), 40)
            temp_tol2 = m_sumMMMMrs(tempN, q1, q2, lambda, muL, muR, vmode, epsilond, T, fc)
            rd = temp_tol == 0.0 && temp_tol2 == 0.0 ? 0.0 :
                 temp_tol == 0.0 ? Inf :
                 abs(log(complex(temp_tol2 / temp_tol)))
        end

        result[qq2] = temp_tol2
    end

    return result
end

"""
    sumMMMMrs11(q1, q2_vec, lambda, muL, muR, vmode, epsilond, T, fc)

Converged double cotunneling sum for n=1→1.
Matches MATLAB `sumMMMMrs11.m`.
"""
function sumMMMMrs11(q1::Int, q2_vec::Vector{Int}, lambda::Float64,
                     muL::Float64, muR::Float64, vmode::Float64, epsilond::Float64,
                     T::Float64, fc::FCCache)
    nq2 = length(q2_vec)
    result = zeros(nq2)

    for qq2 in 1:nq2
        q2 = q2_vec[qq2]
        tempN = round(Int, lambda^2 * 4)
        eps_conv = 1e-14

        temp_tol = m_sumMMMMrs11(tempN, q1, q2, lambda, muL, muR, vmode, epsilond, T, fc)
        tempN += 5
        # MATLAB: m_sumMMMMrs11(max(tempN-1,2), ...)
        temp_tol2 = m_sumMMMMrs11(max(tempN - 1, 2), q1, q2, lambda, muL, muR, vmode, epsilond, T, fc)

        rd = temp_tol == 0.0 && temp_tol2 == 0.0 ? 0.0 :
             temp_tol == 0.0 ? Inf :
             abs(log(complex(temp_tol2 / temp_tol)))

        while isfinite(rd) && rd > eps_conv
            tempN += min(max(round(Int, tempN * 0.5), 20), 40)
            temp_tol = temp_tol2
            temp_tol2 = m_sumMMMMrs11(tempN, q1, q2, lambda, muL, muR, vmode, epsilond, T, fc)
            rd = temp_tol == 0.0 && temp_tol2 == 0.0 ? 0.0 :
                 temp_tol == 0.0 ? Inf :
                 abs(log(complex(temp_tol2 / temp_tol)))
        end

        result[qq2] = temp_tol2
    end

    return result
end
