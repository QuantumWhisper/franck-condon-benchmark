"""Analytically regularized cotunneling integrals I and J.

Uses scipy.special.digamma for complex digamma (fast C code, handles complex natively).
Implements trigamma via asymptotic series with Numba JIT (since scipy.special.polygamma
may not support complex arguments in all versions).

The asymptotic series algorithm matches the Julia port's approach:
- Reflection formula for Re(z) <= 0
- Recurrence shift until |z| >= 20
- Asymptotic expansion with 20 even Bernoulli numbers

Matches MATLAB regularizedI.m, regularizedJ.m, and digammaFcn.m.
"""

import math
import cmath
import numpy as np
import numba
from scipy.special import digamma as scipy_digamma
from .constants import KB_EV
from .fermi_bose import bose_fcn


# ============================================================================
# Bernoulli numbers B_{2k} for k=1..20 (precomputed as float64)
# Used in the asymptotic expansion of trigamma.
# ============================================================================
_BERNOULLI_EVEN = np.array([
    1.0 / 6.0,
    -1.0 / 30.0,
    1.0 / 42.0,
    -1.0 / 30.0,
    5.0 / 66.0,
    -691.0 / 2730.0,
    7.0 / 6.0,
    -3617.0 / 510.0,
    43867.0 / 798.0,
    -174611.0 / 330.0,
    854513.0 / 138.0,
    -236364091.0 / 2730.0,
    8553103.0 / 6.0,
    -23749461029.0 / 870.0,
    8615841276005.0 / 14322.0,
    -7709321041217.0 / 510.0,
    2577687858367.0 / 6.0,
    -26315271553053477373.0 / 1919190.0,
    2929993913841559.0 / 6.0,
    -261082718496449122051.0 / 13530.0,
], dtype=np.float64)

_N_BERNOULLI = 20


# ============================================================================
# Trigamma function via asymptotic series (Numba JIT'd)
# ============================================================================

@numba.njit(cache=True)
def _trigamma_scalar(z):
    """Trigamma function psi'(z) for a single complex128 value.

    Uses recurrence shift + asymptotic series with Bernoulli numbers.
    Handles reflection for Re(z) <= 0.

    psi'(z) = 1/z + 1/(2z^2) + sum_{k=1}^{K} B_{2k} / z^{2k+1}
    """
    z_work = z
    reflection = False

    # Reflection formula for Re(z) <= 0:
    # psi'(z) = (pi/sin(pi*z))^2 - psi'(1-z)
    if z_work.real <= 0:
        reflection = True
        z_work = 1.0 - z_work

    # Recurrence shift: psi'(z) = psi'(z+1) + 1/z^2 until |z| >= 20
    result = 0.0 + 0.0j
    while abs(z_work) < 20:
        result += 1.0 / (z_work * z_work)
        z_work += 1.0

    # Asymptotic expansion
    result += 1.0 / z_work + 1.0 / (2.0 * z_work * z_work)
    z_sq = z_work * z_work
    power = z_sq * z_work  # z^3
    for k in range(_N_BERNOULLI):
        result += _BERNOULLI_EVEN[k] / power
        power *= z_sq

    if reflection:
        sinval = cmath.sin(math.pi * z)
        return (math.pi / sinval) ** 2 - result
    else:
        return result


@numba.njit(cache=True)
def _trigamma_array_flat(z_flat, out):
    """Apply trigamma element-wise to a flat complex array."""
    for i in range(z_flat.shape[0]):
        out[i] = _trigamma_scalar(z_flat[i])


def trigamma(z_arr):
    """Trigamma function for complex array of any shape.

    Parameters
    ----------
    z_arr : ndarray of complex128
        Input array.

    Returns
    -------
    ndarray of complex128
        Trigamma values, same shape as input.
    """
    shape = z_arr.shape
    z_flat = np.ascontiguousarray(z_arr.ravel(), dtype=np.complex128)
    out = np.empty_like(z_flat)
    _trigamma_array_flat(z_flat, out)
    return out.reshape(shape)


# ============================================================================
# Regularized cotunneling integrals
# ============================================================================

def regularized_I(E1, E2, epsilon1_in, epsilon2_in, T):
    """Analytically regularized cotunneling integral I.

    Matches MATLAB regularizedI.m exactly, including broadcasting behavior.
    epsilon1 becomes column, epsilon2 becomes row.

    Parameters
    ----------
    E1 : float
        First energy (typically muL).
    E2 : float
        Second energy (typically muR - (q1-q2)*vmode).
    epsilon1_in : ndarray
        First epsilon array (will be made column).
    epsilon2_in : ndarray
        Second epsilon array (will be made row).
    T : float
        Temperature in Kelvin.

    Returns
    -------
    ndarray, shape (len(epsilon1), len(epsilon2))
    """
    beta = 1.0 / (KB_EV * T)

    # MATLAB: epsilon1 = epsilon1' (transpose to column for broadcasting)
    eps1 = np.asarray(epsilon1_in, dtype=float).reshape(-1, 1)  # column [N, 1]
    eps2 = np.asarray(epsilon2_in, dtype=float).reshape(1, -1)  # row [1, N]

    a1 = 0.5 + 1j * beta * (E2 - eps1) / (2.0 * np.pi)
    a2 = 0.5 - 1j * beta * (E2 - eps2) / (2.0 * np.pi)
    a3 = 0.5 + 1j * beta * (E1 - eps1) / (2.0 * np.pi)
    a4 = 0.5 - 1j * beta * (E1 - eps2) / (2.0 * np.pi)

    # digamma via scipy (handles complex, vectorized C code)
    t1 = scipy_digamma(a1)
    t2 = scipy_digamma(a2)
    t3 = scipy_digamma(a3)
    t4 = scipy_digamma(a4)

    bose_val = bose_fcn(E2 - E1, T)

    with np.errstate(divide='ignore', invalid='ignore'):
        result = bose_val / (eps1 - eps2) * np.real(t1 - t2 - t3 + t4)
    # NaN occurs on diagonal (eps1==eps2); handled by caller zeroing diagonal
    return result


def regularized_J(E1, E2_in, epsilon_in, T):
    """Analytically regularized cotunneling integral J.

    Matches MATLAB regularizedJ.m exactly, including broadcasting behavior.
    epsilon is transposed internally (column or transposed matrix).
    E2 can be scalar or vector (reshaped to row for broadcasting).

    Parameters
    ----------
    E1 : float
        First energy (typically muL).
    E2_in : float or ndarray
        Second energy (scalar or vector).
    epsilon_in : ndarray
        Epsilon values (vector or matrix).
    T : float
        Temperature in Kelvin.

    Returns
    -------
    ndarray
        J values. Shape depends on input dimensions.
    """
    beta = 1.0 / (KB_EV * T)

    # MATLAB: epsilon = epsilon' (transpose)
    if isinstance(epsilon_in, np.ndarray) and epsilon_in.ndim == 2:
        epsilon = epsilon_in.T.copy()  # transpose matrix
    else:
        epsilon = np.asarray(epsilon_in, dtype=float).reshape(-1, 1)  # vector -> column

    # E2: ensure row shape for broadcasting
    if hasattr(E2_in, '__len__'):
        E2 = np.asarray(E2_in, dtype=float).reshape(1, -1)  # row
    else:
        E2 = E2_in  # scalar

    a1 = 0.5 + 1j * beta * (E2 - epsilon) / (2.0 * np.pi)
    a2 = 0.5 + 1j * beta * (E1 - epsilon) / (2.0 * np.pi)

    # trigamma via our Numba implementation (handles complex arrays)
    t1 = trigamma(np.asarray(a1, dtype=np.complex128))
    t2 = trigamma(np.asarray(a2, dtype=np.complex128))

    bose_val = bose_fcn(np.atleast_1d(E2) - E1, T)

    with np.errstate(invalid='ignore'):
        result = beta / (2.0 * np.pi) * bose_val * np.imag(t1 - t2)
    return result
