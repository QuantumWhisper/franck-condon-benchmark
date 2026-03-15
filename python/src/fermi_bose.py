"""Fermi-Dirac and Bose-Einstein distribution functions.

Matches MATLAB fermi.m and BoseFcn.m exactly.
"""

import numpy as np
from .constants import KB_EV


def fermi(x, mu, T):
    """Fermi-Dirac distribution function.

    Parameters
    ----------
    x : float or ndarray
        Energy.
    mu : float
        Chemical potential.
    T : float
        Temperature in Kelvin.

    Returns
    -------
    float or ndarray
        f(x) = 1 / (exp((x - mu) / (kB*T)) + 1).
    """
    with np.errstate(over='ignore'):
        return 1.0 / (np.exp((x - mu) / (KB_EV * T)) + 1.0)


def bose_fcn(x, T):
    """Bose-Einstein distribution function.

    Parameters
    ----------
    x : float or ndarray
        Energy.
    T : float
        Temperature in Kelvin.

    Returns
    -------
    float or ndarray
        n_B(x) = 1 / (exp(x * beta) - 1).
    """
    beta = 1.0 / (KB_EV * T)
    with np.errstate(over='ignore', divide='ignore'):
        return 1.0 / (np.exp(x * beta) - 1.0)
