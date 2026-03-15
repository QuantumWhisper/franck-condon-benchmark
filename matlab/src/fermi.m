function fe = fermi(x,mu,T)
% Fermi-Dirac function

kb = KBoltzmann_ev;
fe = 1./(exp((x-mu)/(kb*T))+1);

end