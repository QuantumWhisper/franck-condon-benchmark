function psi = wavefunction_ho(x,n)
% harmonic oscillator normalized wavefunctions

alpha = 1; % m*\omega/\hbar

psi = 1/sqrt(2^n*factorial(n))*(alpha/pi)^(1/4)*exp(-alpha/2*x.^2).*hermiteH(n,sqrt(alpha)*x);

end