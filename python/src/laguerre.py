"""
Generalized (associated) Laguerre polynomial L_n^alpha(x) via the standard
three-term forward recurrence. Matches MATLAB's laguerreL(n, alpha, x).

    L_0^alpha(x) = 1
    L_1^alpha(x) = 1 + alpha - x
    L_{k+1}^alpha(x) = ((2k + 1 + alpha - x) * L_k - (k + alpha) * L_{k-1}) / (k + 1)
"""

import numba


@numba.njit(cache=True)
def laguerre_L(n, alpha, x):
    """Generalized Laguerre polynomial L_n^alpha(x) via three-term recurrence.

    Parameters
    ----------
    n : int
        Polynomial degree (>= 0).
    alpha : int
        Generalized order.
    x : float
        Evaluation point.

    Returns
    -------
    float
        Value of L_n^alpha(x).
    """
    if n == 0:
        return 1.0
    L_prev = 1.0
    L_curr = 1.0 + alpha - x
    if n == 1:
        return L_curr
    for k in range(1, n):
        L_next = ((2 * k + 1 + alpha - x) * L_curr - (k + alpha) * L_prev) / (k + 1)
        L_prev = L_curr
        L_curr = L_next
    return L_curr
