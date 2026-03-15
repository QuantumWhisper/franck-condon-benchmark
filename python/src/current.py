"""Current computation from steady-state occupation probabilities.

Computes sequential and cotunneling currents from the rate equation solution.
Matches MATLAB current_from_rate_equations.m exactly.

CRITICAL sign conventions:
- n=0->1 (empty to occupied): diffW = w_R - w_L
- n=1->0 (occupied to empty): diffW = w_L - w_R (OPPOSITE!)
- cotunneling: diffW = w_RL - w_LR
"""

import numpy as np
from .matrix import generate_matrix_W
from .solver import solve_steady_state
from .rate import rateW_from_store


def psub_(P, n, q, N):
    """Extract P^n_q from the state vector P.

    Matches MATLAB psub_ in current_from_rate_equations.m.

    Parameters
    ----------
    P : ndarray, shape (2N,)
        State vector.
    n : int
        Charge state (0 or 1).
    q : int
        Phonon state (0-based).
    N : int
        Number of phonon states.

    Returns
    -------
    float
        Occupation probability P^n_q.
    """
    base = 0 if n == 0 else N
    return P[q + base]


def current_from_rate_equations(Vsd, N, vmode, alphaL, alphaR, lambda_, T, eta, Vg, tau, store):
    """Compute sequential and cotunneling currents from steady-state probabilities.

    Matches MATLAB current_from_rate_equations.m exactly.

    Parameters
    ----------
    Vsd : float
        Source-drain bias voltage.
    N : int
        Number of phonon states.
    vmode, alphaL, alphaR, lambda_ : float
        Physical parameters.
    T : float
        Temperature in Kelvin.
    eta : float
        Voltage division factor.
    Vg : float
        Gate voltage.
    tau : float
        Phonon relaxation time.
    store : RateStore
        Pre-computed rates.

    Returns
    -------
    tuple of (float, float, float)
        (I_tol, I_seq, I_cot) — total, sequential, and cotunneling currents.
    """
    # Build and solve rate equation matrix
    M = generate_matrix_W(store, N, vmode, T, tau)
    P = solve_steady_state(M)

    leadR = -1
    leadL = 1
    q2_vec = list(range(N))

    # --- Sequential current for n=0->1 ---
    I_seq0 = 0.0
    for ii in range(N):
        q1 = ii
        w01R = rateW_from_store(store, 0, 1, q1, q2_vec, leadR)
        w01L = rateW_from_store(store, 0, 1, q1, q2_vec, leadL)
        diffW01 = w01R - w01L  # CRITICAL: R - L for n=0->1
        p = psub_(P, 0, q1, N)
        I_seq0 += np.sum(p * diffW01)

    # --- Sequential current for n=1->0 ---
    I_seq1 = 0.0
    for ii in range(N):
        q1 = ii
        w10R = rateW_from_store(store, 1, 0, q1, q2_vec, leadR)
        w10L = rateW_from_store(store, 1, 0, q1, q2_vec, leadL)
        diffW01 = w10L - w10R  # CRITICAL: L - R for n=1->0 (OPPOSITE sign)
        p = psub_(P, 1, q1, N)
        I_seq1 += np.sum(p * diffW01)

    I_seq = I_seq0 + I_seq1

    # --- Cotunneling current ---
    I_cot = 0.0
    for nn in range(2):
        for q1 in range(N):
            p_nn_q1 = psub_(P, nn, q1, N)
            wRL = rateW_from_store(store, nn, nn, q1, q2_vec, leadR)
            wLR = rateW_from_store(store, nn, nn, q1, q2_vec, leadL)
            diffW_nn = wRL - wLR  # cotunneling: RL - LR
            I_cot += np.sum(p_nn_q1 * diffW_nn)

    I_tol = I_seq + I_cot
    return (I_tol, I_seq, I_cot)
