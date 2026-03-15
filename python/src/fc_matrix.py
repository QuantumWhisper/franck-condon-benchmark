"""Franck-Condon matrix elements via generalized Laguerre polynomials.

Matches MATLAB FCMatrixSingle.m and FCMatrix.m exactly.
The FC element <q2|D(lambda)|q1> is computed using:
    sign(q2-q1)^(q1-q2) * lambda^(Q-q) * exp(-lambda^2/2) * sqrt(q!/Q!) * L_q^{Q-q}(lambda^2)
where q = min(q1,q2), Q = max(q1,q2).
"""

import math
import numpy as np
import numba
from .laguerre import laguerre_L


@numba.njit(cache=True)
def fc_matrix_single(q1, q2, lam):
    """Single Franck-Condon matrix element <q2|D(lambda)|q1>.

    Uses log-space computation to avoid factorial overflow for large q values.
    Matches MATLAB FCMatrixSingle.m exactly.

    Parameters
    ----------
    q1, q2 : int
        Phonon state indices.
    lam : float
        Electron-phonon coupling strength lambda.

    Returns
    -------
    float
        FC matrix element M(q1, q2, lambda).
    """
    q = min(q1, q2)
    Q = max(q1, q2)

    l = laguerre_L(q, Q - q, lam * lam)

    # sign(q2-q1)^(q1-q2):
    #   q1==q2: 0^0 = 1
    #   q2>q1:  1^(negative) = 1
    #   q2<q1:  (-1)^(positive integer)
    if q1 == q2:
        sf = 1.0
    elif q2 > q1:
        sf = 1.0
    else:  # q2 < q1
        sf = 1.0 if (q1 - q2) % 2 == 0 else -1.0

    # Handle lambda = 0: FC(q1, q2) = delta(q1, q2)
    if lam == 0.0:
        if Q == q:
            return l  # L_q^0(0) * 1 * 1 * 1
        else:
            return 0.0

    # Use log-space to avoid overflow from factorial and lambda^n
    # M = sf * lambda^(Q-q) * exp(-lambda^2/2) * sqrt(q!/Q!) * L
    # log(|M/L|) = (Q-q)*log(lambda) - lambda^2/2 + 0.5*(lgamma(q+1) - lgamma(Q+1))
    log_coeff = ((Q - q) * math.log(lam)
                 - lam * lam / 2.0
                 + 0.5 * (math.lgamma(q + 1) - math.lgamma(Q + 1)))

    M = sf * math.exp(log_coeff) * l

    if math.isnan(M) or math.isinf(M):
        return 0.0
    return M


class FCCache:
    """FC matrix cache for a fixed lambda value.

    Avoids recomputation of FC elements. Dict-based cache keyed by (q1, q2).
    Replaces MATLAB's memoize(@FCMatrixSingle).
    """

    def __init__(self, lam):
        self._cache = {}
        self.lam = lam

    def get(self, q1, q2):
        """Get FC element for (q1, q2), computing if not cached."""
        key = (q1, q2)
        if key not in self._cache:
            self._cache[key] = fc_matrix_single(q1, q2, self.lam)
        return self._cache[key]

    def matrix(self, q1_range, q2_range):
        """Compute FC matrix M[i,j] = FC(q1_range[i], q2_range[j], lambda).

        Matches MATLAB FCMatrix(q1, q2, lambda).

        Parameters
        ----------
        q1_range : list of int
            Row phonon indices.
        q2_range : list of int
            Column phonon indices.

        Returns
        -------
        ndarray, shape (len(q1_range), len(q2_range))
        """
        nq1 = len(q1_range)
        nq2 = len(q2_range)
        M = np.empty((nq1, nq2))
        for i in range(nq1):
            for j in range(nq2):
                M[i, j] = self.get(q1_range[i], q2_range[j])
        return M
