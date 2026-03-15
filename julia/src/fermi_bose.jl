"""
    fermi(x, mu, T)

Fermi-Dirac distribution function. Matches MATLAB `fermi.m`.
"""
function fermi(x, mu, T)
    return 1.0 / (exp((x - mu) / (KB_EV * T)) + 1.0)
end

"""
    bose_fcn(x, T)

Bose-Einstein distribution function. Matches MATLAB `BoseFcn.m`.
"""
function bose_fcn(x, T)
    beta = 1.0 / (KB_EV * T)
    return 1.0 / (exp(x * beta) - 1.0)
end
