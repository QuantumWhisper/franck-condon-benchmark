"""Cotunneling sum functions with adaptive convergence wrappers.

Inner functions (m_sumMMr, m_sumMMMMrs, etc.) compute the cotunneling sums for
a given truncation N. Outer wrappers adaptively increase N until convergence.

Matches MATLAB sumMMr.m, sumMMMMrs.m, sumMMr11.m, sumMMMMrs11.m and their
m_ inner functions exactly, including convergence criteria:
- sumMMr, sumMMr11: log10 convergence criterion
- sumMMMMrs, sumMMMMrs11: log (natural) convergence criterion
"""

import numpy as np
from .regularized import regularized_I, regularized_J


# ===================== Inner computation functions =====================

def m_sumMMr(N, q1, q2_vec, lambda_, muL, muR, vmode, epsilond, T, fc):
    """Single cotunneling sum over virtual states r for n=0->0 transition.

    Matches MATLAB m_sumMMr.m exactly.

    Parameters
    ----------
    N : int
        Truncation for the sum over r.
    q1 : int
        Initial phonon state.
    q2_vec : list of int
        Final phonon states.
    lambda_ : float
        Electron-phonon coupling.
    muL, muR : float
        Chemical potentials of left and right leads.
    vmode : float
        Phonon energy in eV.
    epsilond : float
        Molecular level energy.
    T : float
        Temperature in Kelvin.
    fc : FCCache
        Franck-Condon matrix cache.

    Returns
    -------
    ndarray, shape (nq2,)
        Sum values for each q2.
    """
    rr = list(range(N))
    nq2 = len(q2_vec)

    # FCMatrix(q2, rr)' .* FCMatrix(q1, rr)' -> [nr, nq2]
    fc_q2_rr = fc.matrix(q2_vec, rr)    # [nq2, nr]
    fc_q1_rr = fc.matrix([q1], rr)       # [1, nr]

    # MATLAB: abs(FCMatrix(q2,rr,lambda)'.*FCMatrix(q1,rr,lambda)').^2
    MM_square = np.abs(fc_q2_rr.T * fc_q1_rr.T) ** 2  # [nr, nq2]

    # Jr = regularizedJ(muL, muR-(q1-q2)*vmode, epsilond-(q1-rr)*vmode, T)
    E2_vec = muR - (q1 - np.array(q2_vec)) * vmode   # [nq2]
    eps_vec = epsilond - (q1 - np.array(rr)) * vmode  # [nr]
    Jr = regularized_J(muL, E2_vec, eps_vec, T)        # [nr, nq2]

    with np.errstate(invalid='ignore'):
        target = MM_square * Jr
    target = np.where(np.isfinite(target), target, 0.0)
    return np.sum(target, axis=0)  # [nq2]


def m_sumMMr11(N, q1, q2_vec, lambda_, muL, muR, vmode, epsilond, T, fc):
    """Single cotunneling sum for n=1->1 transition.

    Matches MATLAB m_sumMMr11.m exactly.
    Key difference from m_sumMMr: epsilon uses (q2-rr) instead of (q1-rr).
    """
    rr = list(range(N))
    nq2 = len(q2_vec)

    # MATLAB: abs(FCMatrix(q2,rr).*FCMatrix(q1,rr)).^2'
    fc_q2_rr = fc.matrix(q2_vec, rr)    # [nq2, nr]
    fc_q1_rr = fc.matrix([q1], rr)       # [1, nr]
    MM_square = (np.abs(fc_q2_rr * fc_q1_rr) ** 2).T  # [nr, nq2]

    # epsilon = epsilond + (q2' - rr) * vmode -> matrix [nq2, nr]
    q2_col = np.array(q2_vec).reshape(-1, 1)
    rr_row = np.array(rr).reshape(1, -1)
    eps_matrix = epsilond + (q2_col - rr_row) * vmode  # [nq2, nr]

    E2_vec = muR - (q1 - np.array(q2_vec)) * vmode
    Jr = regularized_J(muL, E2_vec, eps_matrix, T)     # [nr, nq2] after internal transpose

    with np.errstate(invalid='ignore'):
        target = MM_square * Jr
    target = np.where(np.isfinite(target), target, 0.0)
    return np.sum(target, axis=0)  # [nq2]


def m_sumMMMMrs(N, q1, q2, lambda_, muL, muR, vmode, epsilond, T, fc):
    """Double cotunneling sum over r,s for n=0->0 transition.

    Matches MATLAB m_sumMMMMrs.m + m_MMMM_quad_Irs.m.
    Returns scalar.
    """
    rr = list(range(N))
    ss = list(range(N))
    nr = N

    # MMMM_quad2D = FC(q2,rr)' .* conj(FC(q1,rr))' .* conj(FC(q2,ss)) .* FC(q1,ss)
    fc_q2_rr = fc.matrix([q2], rr)  # [1, nr]
    fc_q1_rr = fc.matrix([q1], rr)  # [1, nr]
    fc_q2_ss = fc.matrix([q2], ss)  # [1, nr]
    fc_q1_ss = fc.matrix([q1], ss)  # [1, nr]

    # FC values are real, so conj is identity. But keep conj for correctness.
    MMMM = (fc_q2_rr.T * np.conj(fc_q1_rr.T)
            * np.conj(fc_q2_ss) * fc_q1_ss)  # [nr, nr]

    # Irs = regularizedI(muL, muR-(q1-q2)*vmode, epsilond-(q1-rr)*vmode, epsilond-(q1-ss)*vmode, T)
    E2 = muR - (q1 - q2) * vmode
    eps1 = epsilond - (q1 - np.array(rr)) * vmode
    eps2 = epsilond - (q1 - np.array(ss)) * vmode
    Irs = regularized_I(muL, E2, eps1, eps2, T)  # [nr, nr]

    with np.errstate(invalid='ignore'):
        target = MMMM * Irs
    # Zero diagonal
    np.fill_diagonal(target, 0.0)
    target = np.where(np.isfinite(target), target, 0.0)
    return np.sum(target)


def m_sumMMMMrs11(N, q1, q2, lambda_, muL, muR, vmode, epsilond, T, fc):
    """Double cotunneling sum for n=1->1 transition.

    Matches MATLAB m_sumMMMMrs11.m.
    Key difference from n=0->0: epsilon uses epsilond+(q2-rr)*vmode.
    """
    rr = list(range(N))
    ss = list(range(N))
    nr = N

    fc_q2_rr = fc.matrix([q2], rr)
    fc_q1_rr = fc.matrix([q1], rr)
    fc_q2_ss = fc.matrix([q2], ss)
    fc_q1_ss = fc.matrix([q1], ss)

    MMMM = (fc_q2_rr.T * np.conj(fc_q1_rr.T)
            * np.conj(fc_q2_ss) * fc_q1_ss)

    # Different epsilon for n=1->1: epsilond+(q2-rr)*vmode
    E2 = muR - (q1 - q2) * vmode
    eps1 = epsilond + (q2 - np.array(rr)) * vmode
    eps2 = epsilond + (q2 - np.array(ss)) * vmode
    Irs = regularized_I(muL, E2, eps1, eps2, T)

    with np.errstate(invalid='ignore'):
        target = MMMM * Irs
    np.fill_diagonal(target, 0.0)
    target = np.where(np.isfinite(target), target, 0.0)
    return np.sum(target)


# ===================== Convergence wrappers =====================

def sumMMr(q1, q2_vec, lambda_, muL, muR, vmode, epsilond, T, fc):
    """Converged single cotunneling sum for n=0->0.

    Uses adaptive N with log10 convergence criterion.
    Matches MATLAB sumMMr.m.
    """
    tempN = round(lambda_ ** 2.2 * 3)
    epsilon = 1e-14

    temp_tol = m_sumMMr(tempN, q1, q2_vec, lambda_, muL, muR, vmode, epsilond, T, fc)
    tempN += 5
    temp_tol2 = m_sumMMr(tempN, q1, q2_vec, lambda_, muL, muR, vmode, epsilond, T, fc)

    # MATLAB uses log10 for convergence (not log)
    with np.errstate(divide='ignore', invalid='ignore'):
        relative_diff = np.abs(np.log10(temp_tol2 / temp_tol + 0j).real)
        relative_diff = np.where(np.isfinite(relative_diff), relative_diff, 0.0)

    while np.any(relative_diff > epsilon):
        locs = relative_diff > epsilon
        tempN += min(max(round(tempN * 0.5), 10), 20)
        temp_tol = temp_tol2.copy()
        q2_subset = [q2_vec[i] for i in range(len(q2_vec)) if locs[i]]
        temp_tol2_partial = m_sumMMr(tempN, q1, q2_subset, lambda_, muL, muR, vmode, epsilond, T, fc)
        temp_tol2[locs] = temp_tol2_partial
        with np.errstate(divide='ignore', invalid='ignore'):
            relative_diff = np.abs(np.log10(temp_tol2 / temp_tol + 0j).real)
            relative_diff = np.where(np.isfinite(relative_diff), relative_diff, 0.0)

    return temp_tol2


def sumMMr11(q1, q2_vec, lambda_, muL, muR, vmode, epsilond, T, fc):
    """Converged single cotunneling sum for n=1->1.

    Matches MATLAB sumMMr11.m.
    """
    tempN = round(lambda_ ** 2.2 * 3)
    epsilon = 1e-14

    temp_tol = m_sumMMr11(tempN, q1, q2_vec, lambda_, muL, muR, vmode, epsilond, T, fc)
    tempN += 5
    temp_tol2 = m_sumMMr11(tempN, q1, q2_vec, lambda_, muL, muR, vmode, epsilond, T, fc)

    with np.errstate(divide='ignore', invalid='ignore'):
        relative_diff = np.abs(np.log10(temp_tol2 / temp_tol + 0j).real)
        relative_diff = np.where(np.isfinite(relative_diff), relative_diff, 0.0)

    while np.any(relative_diff > epsilon):
        locs = relative_diff > epsilon
        tempN += min(max(round(tempN * 0.5), 10), 20)
        temp_tol = temp_tol2.copy()
        q2_subset = [q2_vec[i] for i in range(len(q2_vec)) if locs[i]]
        temp_tol2_partial = m_sumMMr11(tempN, q1, q2_subset, lambda_, muL, muR, vmode, epsilond, T, fc)
        temp_tol2[locs] = temp_tol2_partial
        with np.errstate(divide='ignore', invalid='ignore'):
            relative_diff = np.abs(np.log10(temp_tol2 / temp_tol + 0j).real)
            relative_diff = np.where(np.isfinite(relative_diff), relative_diff, 0.0)

    return temp_tol2


def sumMMMMrs(q1, q2_vec, lambda_, muL, muR, vmode, epsilond, T, fc):
    """Converged double cotunneling sum for n=0->0.

    Processes each q2 individually. Uses log (natural) convergence criterion.
    Matches MATLAB sumMMMMrs.m.
    """
    nq2 = len(q2_vec)
    result = np.zeros(nq2)

    for qq2 in range(nq2):
        q2 = q2_vec[qq2]
        tempN = round(lambda_ ** 2 * 4)
        eps_conv = 1e-14

        temp_tol = m_sumMMMMrs(tempN, q1, q2, lambda_, muL, muR, vmode, epsilond, T, fc)
        tempN += 5
        temp_tol2 = m_sumMMMMrs(tempN, q1, q2, lambda_, muL, muR, vmode, epsilond, T, fc)

        # MATLAB uses log (natural) for convergence
        if temp_tol == 0.0 and temp_tol2 == 0.0:
            rd = 0.0
        elif temp_tol == 0.0:
            rd = float('inf')
        else:
            rd = abs(np.log(complex(temp_tol2 / temp_tol)))

        while np.isfinite(rd) and rd > eps_conv:
            temp_tol = temp_tol2
            tempN += min(max(round(tempN * 0.5), 20), 40)
            temp_tol2 = m_sumMMMMrs(tempN, q1, q2, lambda_, muL, muR, vmode, epsilond, T, fc)
            if temp_tol == 0.0 and temp_tol2 == 0.0:
                rd = 0.0
            elif temp_tol == 0.0:
                rd = float('inf')
            else:
                rd = abs(np.log(complex(temp_tol2 / temp_tol)))

        result[qq2] = temp_tol2

    return result


def sumMMMMrs11(q1, q2_vec, lambda_, muL, muR, vmode, epsilond, T, fc):
    """Converged double cotunneling sum for n=1->1.

    Matches MATLAB sumMMMMrs11.m.
    """
    nq2 = len(q2_vec)
    result = np.zeros(nq2)

    for qq2 in range(nq2):
        q2 = q2_vec[qq2]
        tempN = round(lambda_ ** 2 * 4)
        eps_conv = 1e-14

        temp_tol = m_sumMMMMrs11(tempN, q1, q2, lambda_, muL, muR, vmode, epsilond, T, fc)
        tempN += 5
        # MATLAB: m_sumMMMMrs11(max(tempN-1,2), ...)
        temp_tol2 = m_sumMMMMrs11(max(tempN - 1, 2), q1, q2, lambda_, muL, muR, vmode, epsilond, T, fc)

        if temp_tol == 0.0 and temp_tol2 == 0.0:
            rd = 0.0
        elif temp_tol == 0.0:
            rd = float('inf')
        else:
            rd = abs(np.log(complex(temp_tol2 / temp_tol)))

        while np.isfinite(rd) and rd > eps_conv:
            tempN += min(max(round(tempN * 0.5), 20), 40)
            temp_tol = temp_tol2
            temp_tol2 = m_sumMMMMrs11(tempN, q1, q2, lambda_, muL, muR, vmode, epsilond, T, fc)
            if temp_tol == 0.0 and temp_tol2 == 0.0:
                rd = 0.0
            elif temp_tol == 0.0:
                rd = float('inf')
            else:
                rd = abs(np.log(complex(temp_tol2 / temp_tol)))

        result[qq2] = temp_tol2

    return result
