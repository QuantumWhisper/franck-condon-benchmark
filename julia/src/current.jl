"""
    psub_(P, n, q, N)

Extract P^n_q from the state vector P.
Matches MATLAB `psub_` in current_from_rate_equations.m.
"""
@inline function psub_(P::Vector{Float64}, n::Int, q::Int, N::Int)
    base = n == 0 ? 0 : N
    @inbounds return P[q + 1 + base]
end

"""
    current_from_rate_equations(Vsd, N, vmode, alphaL, alphaR, lambda, T, eta, Vg, tau, store)

Compute sequential and cotunneling currents from steady-state probabilities.
Matches MATLAB `current_from_rate_equations.m` exactly.

Returns (I_tol, I_seq, I_cot).
"""
function current_from_rate_equations(Vsd::Float64, N::Int, vmode::Float64,
                                      alphaL::Float64, alphaR::Float64, lambda::Float64,
                                      T::Float64, eta::Float64, Vg::Float64, tau::Float64,
                                      store::RateStore)
    # Build and solve rate equation matrix
    M = generate_matrix_W(store, N, vmode, T, tau)
    P = solve_steady_state(M)

    leadR = -1
    leadL = 1
    q2_vec = collect(0:N-1)

    # --- Sequential current for n=0→1 ---
    I_seq0 = 0.0
    for ii in 1:N
        q1 = ii - 1
        w01R = rateW_from_store(store, 0, 1, q1, q2_vec, leadR)
        w01L = rateW_from_store(store, 0, 1, q1, q2_vec, leadL)
        diffW01 = w01R .- w01L
        p = psub_(P, 0, q1, N)
        I_seq0 += sum(p .* diffW01)
    end

    # --- Sequential current for n=1→0 ---
    I_seq1 = 0.0
    for ii in 1:N
        q1 = ii - 1
        w10R = rateW_from_store(store, 1, 0, q1, q2_vec, leadR)
        w10L = rateW_from_store(store, 1, 0, q1, q2_vec, leadL)
        diffW01 = w10L .- w10R  # Note: L - R for n=1→0 (opposite sign)
        p = psub_(P, 1, q1, N)
        I_seq1 += sum(p .* diffW01)
    end

    I_seq = I_seq0 + I_seq1

    # --- Cotunneling current ---
    I_cot = 0.0
    for nn in 0:1
        for q1 in 0:(N-1)
            p_nn_q1 = psub_(P, nn, q1, N)
            wRL = rateW_from_store(store, nn, nn, q1, q2_vec, leadR)
            wLR = rateW_from_store(store, nn, nn, q1, q2_vec, leadL)
            diffW_nn = wRL .- wLR
            I_cot += sum(p_nn_q1 .* diffW_nn)
        end
    end

    I_tol = I_seq + I_cot
    return (I_tol, I_seq, I_cot)
end
