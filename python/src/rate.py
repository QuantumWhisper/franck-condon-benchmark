"""Core rate computation and pre-computation store.

Matches MATLAB m_rateW.m, rateW.m, calculateAllRateW.m exactly.
Rate store replaces MATLAB's memoize cache with a dict keyed by (n1, n2, q1, lead).
"""

import numpy as np
from .constants import HBAR_EV
from .fermi_bose import fermi
from .cotunneling import sumMMr, sumMMMMrs, sumMMr11, sumMMMMrs11


def spin_degeneracy(n1, n2):
    """Spin degeneracy factor s(n1, n2).

    Matches MATLAB s_ in m_rateW.m.
    - s=2 for n1=0 (empty to occupied, spin up or down)
    - s=1 for n1=1, n2=0 (occupied to empty)
    - s=2 for n1=1, n2=1 (cotunneling)
    """
    if n1 == 0:
        return 2.0
    else:  # n1 == 1
        return 1.0 if n2 == 0 else 2.0


def m_rateW(n1, n2, q1, q2_vec, vmode, alphaL, alphaR, lambda_, Vsd, T, eta, lead, Vg, fc):
    """Core rate computation. Matches MATLAB m_rateW.m exactly.

    Returns vector of rates for each q2 value.

    Parameters
    ----------
    n1, n2 : int
        Initial and final charge states (0 or 1).
    q1 : int
        Initial phonon state.
    q2_vec : list of int
        Final phonon states.
    vmode : float
        Phonon energy in eV.
    alphaL, alphaR : float
        Tunnel coupling ratios.
    lambda_ : float
        Electron-phonon coupling.
    Vsd : float
        Source-drain bias voltage.
    T : float
        Temperature in Kelvin.
    eta : float
        Voltage division factor.
    lead : int
        Lead index (+1 for left, -1 for right).
    Vg : float
        Gate voltage.
    fc : FCCache
        Franck-Condon matrix cache.

    Returns
    -------
    ndarray, shape (len(q2_vec),)
        Transition rates.
    """
    gammaL = alphaL * vmode
    gammaR = alphaR * vmode
    epsilond = 0.0 + Vg
    muL = eta * Vsd
    muR = -(1.0 - eta) * Vsd

    s = spin_degeneracy(n1, n2)
    gamma = gammaL if lead == 1 else gammaR
    mu = muL if lead == 1 else muR

    q2_arr = np.array(q2_vec, dtype=float)

    if n1 == 1 and n2 == 0:
        # Sequential tunneling: 1->0 (occupied to empty)
        fc_vals = np.array([fc.get(q1, q) for q in q2_vec])
        w = (s * gamma / HBAR_EV * fc_vals ** 2
             * (1.0 - fermi(epsilond - (q2_arr - q1) * vmode, mu, T)))

    elif n1 == 0 and n2 == 1:
        # Sequential tunneling: 0->1 (empty to occupied)
        fc_vals = np.array([fc.get(q1, q) for q in q2_vec])
        w = (s * gamma / HBAR_EV * fc_vals ** 2
             * fermi(epsilond + (q2_arr - q1) * vmode, mu, T))

    elif n1 == 0 and n2 == 0:
        # Cotunneling: 0->0
        if lead == 1:  # L -> R
            w = (s / (2.0 * np.pi * HBAR_EV) * gammaL * gammaR
                 * (sumMMr(q1, q2_vec, lambda_, muL, muR, vmode, epsilond, T, fc)
                    + sumMMMMrs(q1, q2_vec, lambda_, muL, muR, vmode, epsilond, T, fc)))
        else:  # R -> L
            w = (s / (2.0 * np.pi * HBAR_EV) * gammaL * gammaR
                 * (sumMMr(q1, q2_vec, lambda_, muR, muL, vmode, epsilond, T, fc)
                    + sumMMMMrs(q1, q2_vec, lambda_, muR, muL, vmode, epsilond, T, fc)))

    else:  # n1 == 1 and n2 == 1
        # Cotunneling: 1->1
        if lead == 1:
            w = (s / (2.0 * np.pi * HBAR_EV) * gammaL * gammaR
                 * (sumMMr11(q1, q2_vec, lambda_, muL, muR, vmode, epsilond, T, fc)
                    + sumMMMMrs11(q1, q2_vec, lambda_, muL, muR, vmode, epsilond, T, fc)))
        else:
            w = (s / (2.0 * np.pi * HBAR_EV) * gammaL * gammaR
                 * (sumMMr11(q1, q2_vec, lambda_, muR, muL, vmode, epsilond, T, fc)
                    + sumMMMMrs11(q1, q2_vec, lambda_, muR, muL, vmode, epsilond, T, fc)))

    # Replace NaN with 0 (matches MATLAB: w(isnan(w)) = 0)
    w = np.where(np.isnan(w), 0.0, w)
    return w


class RateStore(dict):
    """Pre-computed rates indexed by (n1, n2, q1, lead) -> ndarray.

    Replaces MATLAB's memoize cache. Values are numpy arrays of length N
    containing rates for all q2 values.
    """
    pass


def calculate_all_rateW(store, N, vmode, alphaL, alphaR, lambda_, Vsd, T, eta, lead, Vg, fc):
    """Pre-compute all rates for a given lead and bias point.

    Fills the rate store with rates for all (n1, n2, q1) combinations.
    Matches MATLAB calculateAllRateW.m.

    Parameters
    ----------
    store : RateStore
        Rate store to fill.
    N : int
        Number of phonon states.
    Other parameters match m_rateW.
    """
    q2_vec = list(range(N))
    for n1 in range(2):
        for n2 in range(2):
            for ii in range(N):
                q1 = ii
                w = m_rateW(n1, n2, q1, q2_vec, vmode, alphaL, alphaR, lambda_,
                            Vsd, T, eta, lead, Vg, fc)
                store[(n1, n2, q1, lead)] = w


def rateW_from_store(store, n1, n2, q1, q2_indices, lead):
    """Retrieve pre-computed rates from store.

    Parameters
    ----------
    store : RateStore
    n1, n2, q1, lead : int
        Rate indices.
    q2_indices : int or list of int
        0-based phonon indices to retrieve.

    Returns
    -------
    float or ndarray
    """
    w_full = store[(n1, n2, q1, lead)]
    if isinstance(q2_indices, (int, np.integer)):
        return w_full[q2_indices]
    else:
        return w_full[np.array(q2_indices)]


def rateW_lead_from_store(store, n1, n2, q1, q2):
    """Rate summed over both leads. Matches MATLAB rateW_lead.m.

    Parameters
    ----------
    store : RateStore
    n1, n2, q1 : int
    q2 : int or list of int
        0-based phonon indices.

    Returns
    -------
    float or ndarray
    """
    wl = rateW_from_store(store, n1, n2, q1, q2, 1)
    wr = rateW_from_store(store, n1, n2, q1, q2, -1)
    return wl + wr
