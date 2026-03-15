"""
    peq_(q, vmode, T)

Equilibrium Bose-Einstein phonon occupation probability.
Matches MATLAB `peq_` in generateMatrixW.m.
"""
function peq_(q, vmode::Float64, T::Float64)
    beta = 1.0 / (KB_EV * T)
    return exp.(-q .* vmode .* beta) .* (1.0 .- exp(-vmode * beta))
end

"""
    sigma_W(store, N, n1, n2, q1)

Sum of rates over all q2 for given (n1, n2, q1), both leads combined.
Matches MATLAB sigmaW00, sigmaW01, sigmaW10, sigmaW11 in generateMatrixW.m.
"""
function sigma_W(store::RateStore, N::Int, n1::Int, n2::Int, q1::Int)
    q2_vec = collect(0:N-1)
    w = rateW_lead_from_store(store, n1, n2, q1, q2_vec)
    return sum(w)
end

"""
    generate_matrix_W(store, N, vmode, T, tau)

Build the 2N×2N rate equation matrix W.
Matches MATLAB `generateMatrixW.m` exactly.

State ordering: P^0_0, P^0_1, ..., P^0_{N-1}, P^1_0, P^1_1, ..., P^1_{N-1}
"""
function generate_matrix_W(store::RateStore, N::Int, vmode::Float64, T::Float64, tau::Float64)
    M = zeros(2N, 2N)

    # ---- n = 0 block (rows 1:N) ----
    base = 0
    for ii in 1:N
        m1Index = base + ii
        for jj in 1:(2N)
            modjj = mod(jj, N)
            if modjj == 0
                modjj = N
            end
            P_q_index = modjj - 1  # 0-based phonon index

            if jj <= N
                if ii == modjj
                    # Diagonal of n=0 block
                    n1, n2 = 0, 0
                    q1 = P_q_index
                    q2 = P_q_index
                    M[m1Index, jj] = rateW_lead_from_store(store, n1, n2, q1, q2) -
                                     sigma_W(store, N, 0, 0, q1) -
                                     sigma_W(store, N, 0, 1, q1) -
                                     1.0 / tau + peq_(P_q_index, vmode, T) / tau
                else
                    # Off-diagonal within n=0 block
                    n1, n2 = 0, 0
                    q1 = P_q_index
                    q2 = ii - 1
                    M[m1Index, jj] = rateW_lead_from_store(store, n1, n2, q1, q2) +
                                     peq_(q2, vmode, T) / tau
                end
            else
                # n=1→0 transition block
                n1, n2 = 1, 0
                q1 = P_q_index
                q2 = ii - 1
                M[m1Index, jj] = rateW_lead_from_store(store, n1, n2, q1, q2)
            end
        end
    end

    # ---- n = 1 block (rows N+1:2N) ----
    base = N
    for ii in 1:N
        m1Index = base + ii
        for jj in 1:(2N)
            modjj = mod(jj, N)
            if modjj == 0
                modjj = N
            end

            if jj <= N
                # n=0→1 transition block
                n1, n2 = 0, 1
                q1 = modjj - 1
                q2 = ii - 1
                M[m1Index, jj] = rateW_lead_from_store(store, n1, n2, q1, q2)
            else
                if modjj == ii
                    # Diagonal of n=1 block
                    n1, n2 = 1, 1
                    q1 = modjj - 1
                    q2 = ii - 1
                    M[m1Index, jj] = rateW_lead_from_store(store, n1, n2, q1, q2) -
                                     sigma_W(store, N, 1, 0, q2) -
                                     sigma_W(store, N, 1, 1, q2) -
                                     1.0 / tau + peq_(q2, vmode, T) / tau
                else
                    # Off-diagonal within n=1 block
                    n1, n2 = 1, 1
                    q1 = modjj - 1
                    q2 = ii - 1
                    M[m1Index, jj] = rateW_lead_from_store(store, n1, n2, q1, q2) +
                                     peq_(q2, vmode, T) / tau
                end
            end
        end
    end

    return M
end
