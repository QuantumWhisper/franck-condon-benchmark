"""Main I-V simulation driver.

Matches Julia FranckCondon.jl simulate_iv function and MATLAB run_benchmark.m
simulation loop exactly.
"""

import numpy as np
from .constants import ELEMENTARY_CHARGE
from .fc_matrix import FCCache
from .rate import RateStore, calculate_all_rateW
from .current import current_from_rate_equations


def simulate_iv(N, vmode, alphaL, alphaR, lambda_, Vsd_vec, T, eta, Vg, tau, verbose=True):
    """Run the full I-V simulation.

    Parameters
    ----------
    N : int
        Number of phonon states.
    vmode : float
        Phonon energy in eV.
    alphaL, alphaR : float
        Tunnel coupling ratios (Gamma/hbar*omega).
    lambda_ : float
        Electron-phonon coupling strength.
    Vsd_vec : ndarray
        Bias voltage sweep values.
    T : float
        Temperature in Kelvin.
    eta : float
        Voltage division factor.
    Vg : float
        Gate voltage.
    tau : float
        Phonon relaxation time.
    verbose : bool
        Print progress.

    Returns
    -------
    tuple of (ndarray, ndarray, ndarray, ndarray)
        (Vsd_vec, I_tol, I_seq, I_cot) in SI units (Amperes).
    """
    nVsd = len(Vsd_vec)
    I_tol = np.zeros(nVsd)
    I_seq = np.zeros(nVsd)
    I_cot = np.zeros(nVsd)

    # FC matrix cache persists across all bias points
    fc = FCCache(lambda_)

    for vv in range(nVsd):
        v = Vsd_vec[vv]
        if verbose:
            print(f"Vsd = {v:.4f} V  ({vv + 1}/{nVsd})")

        # Pre-compute rate matrices for both leads
        store = RateStore()
        calculate_all_rateW(store, N, vmode, alphaL, alphaR, lambda_, v, T, eta, 1, Vg, fc)
        calculate_all_rateW(store, N, vmode, alphaL, alphaR, lambda_, v, T, eta, -1, Vg, fc)

        # Compute current
        I_tol[vv], I_seq[vv], I_cot[vv] = current_from_rate_equations(
            v, N, vmode, alphaL, alphaR, lambda_, T, eta, Vg, tau, store)

        # Sign convention: current flows opposite to bias
        if v != 0.0:
            s = -np.sign(v)
            I_tol[vv] *= s
            I_seq[vv] *= s
            I_cot[vv] *= s

    # Convert to SI (Amperes)
    e = ELEMENTARY_CHARGE
    I_tol *= e
    I_seq *= e
    I_cot *= e

    return (Vsd_vec, I_tol, I_seq, I_cot)
