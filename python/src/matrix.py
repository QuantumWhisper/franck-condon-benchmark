"""Rate equation matrix W construction.

Builds the 2N x 2N matrix W for the rate equations 0 = WP.
State ordering: P^0_0, P^0_1, ..., P^0_{N-1}, P^1_0, P^1_1, ..., P^1_{N-1}

Matches MATLAB generateMatrixW.m exactly, including phonon relaxation terms.
"""

import numpy as np
from .constants import KB_EV
from .rate import rateW_lead_from_store


def peq_(q, vmode, T):
    """Equilibrium Bose-Einstein phonon occupation probability.

    Matches MATLAB peq_ in generateMatrixW.m.

    Parameters
    ----------
    q : int or ndarray
        Phonon state index.
    vmode : float
        Phonon energy in eV.
    T : float
        Temperature in Kelvin.

    Returns
    -------
    float or ndarray
    """
    beta = 1.0 / (KB_EV * T)
    return np.exp(-q * vmode * beta) * (1.0 - np.exp(-vmode * beta))


def sigma_W(store, N, n1, n2, q1):
    """Sum of rates over all q2 for given (n1, n2, q1), both leads combined.

    Matches MATLAB sigmaW00, sigmaW01, sigmaW10, sigmaW11 in generateMatrixW.m.
    """
    q2_vec = list(range(N))
    w = rateW_lead_from_store(store, n1, n2, q1, q2_vec)
    return np.sum(w)


def generate_matrix_W(store, N, vmode, T, tau):
    """Build the 2N x 2N rate equation matrix W.

    Matches MATLAB generateMatrixW.m exactly.

    Parameters
    ----------
    store : RateStore
        Pre-computed rates for current bias point.
    N : int
        Number of phonon states.
    vmode : float
        Phonon energy in eV.
    T : float
        Temperature in Kelvin.
    tau : float
        Phonon relaxation time (Inf = unequilibrated).

    Returns
    -------
    ndarray, shape (2N, 2N)
        Rate equation matrix W.
    """
    M = np.zeros((2 * N, 2 * N))

    # ---- n = 0 block (rows 0:N-1) ----
    base = 0
    for ii in range(1, N + 1):
        m1Index = base + ii - 1  # 0-based
        for jj in range(1, 2 * N + 1):
            modjj = jj % N
            if modjj == 0:
                modjj = N
            P_q_index = modjj - 1  # 0-based phonon index

            if jj <= N:
                if ii == modjj:
                    # Diagonal of n=0 block
                    q1 = P_q_index
                    M[m1Index, jj - 1] = (
                        rateW_lead_from_store(store, 0, 0, q1, q1)
                        - sigma_W(store, N, 0, 0, q1)
                        - sigma_W(store, N, 0, 1, q1)
                        - 1.0 / tau + peq_(P_q_index, vmode, T) / tau
                    )
                else:
                    # Off-diagonal within n=0 block
                    q1 = P_q_index
                    q2 = ii - 1
                    M[m1Index, jj - 1] = (
                        rateW_lead_from_store(store, 0, 0, q1, q2)
                        + peq_(q2, vmode, T) / tau
                    )
            else:
                # n=1->0 transition block
                q1 = P_q_index
                q2 = ii - 1
                M[m1Index, jj - 1] = rateW_lead_from_store(store, 1, 0, q1, q2)

    # ---- n = 1 block (rows N:2N-1) ----
    base = N
    for ii in range(1, N + 1):
        m1Index = base + ii - 1  # 0-based
        for jj in range(1, 2 * N + 1):
            modjj = jj % N
            if modjj == 0:
                modjj = N

            if jj <= N:
                # n=0->1 transition block
                q1 = modjj - 1
                q2 = ii - 1
                M[m1Index, jj - 1] = rateW_lead_from_store(store, 0, 1, q1, q2)
            else:
                if modjj == ii:
                    # Diagonal of n=1 block
                    q1 = modjj - 1
                    q2 = ii - 1
                    M[m1Index, jj - 1] = (
                        rateW_lead_from_store(store, 1, 1, q1, q2)
                        - sigma_W(store, N, 1, 0, q2)
                        - sigma_W(store, N, 1, 1, q2)
                        - 1.0 / tau + peq_(q2, vmode, T) / tau
                    )
                else:
                    # Off-diagonal within n=1 block
                    q1 = modjj - 1
                    q2 = ii - 1
                    M[m1Index, jj - 1] = (
                        rateW_lead_from_store(store, 1, 1, q1, q2)
                        + peq_(q2, vmode, T) / tau
                    )

    return M
