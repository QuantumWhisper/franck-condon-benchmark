"""Steady-state solver for the rate equation system 0 = WP.

Matches Julia solver.jl: augmented system approach (replace last row with
normalization constraint sum(P) = 1), solve via least-squares, clamp negatives,
renormalize.
"""

import numpy as np


def solve_steady_state(M):
    """Solve for steady-state occupation probabilities from 0 = WP.

    Uses the augmented system approach:
    1. Replace last row of W with normalization: sum(P) = 1
    2. Solve via numpy.linalg.lstsq
    3. Clamp negatives to 0
    4. Renormalize so sum(P) = 1

    Parameters
    ----------
    M : ndarray, shape (2N, 2N)
        Rate equation matrix W.

    Returns
    -------
    ndarray, shape (2N,)
        Steady-state probability vector P.
    """
    n = M.shape[1]

    # Replace last row with normalization constraint
    M_aug = M.copy()
    M_aug[-1, :] = 1.0
    d = np.zeros(n)
    d[-1] = 1.0

    # Solve via QR (backslash equivalent)
    P, _, _, _ = np.linalg.lstsq(M_aug, d, rcond=None)

    # Clamp negatives and renormalize
    P = np.maximum(P, 0.0)
    s = np.sum(P)
    if s > 0:
        P /= s

    return P
