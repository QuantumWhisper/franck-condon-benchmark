"""
Franck-Condon blockade I-V simulation — Python+Numba implementation.

A faithful 1:1 port of the MATLAB reference, following the Julia port structure.
Reference: Koch, von Oppen, Glazman, PRB 74, 205438 (2006).
"""

from .simulate import simulate_iv
from .constants import ELEMENTARY_CHARGE

__all__ = ["simulate_iv", "ELEMENTARY_CHARGE"]
