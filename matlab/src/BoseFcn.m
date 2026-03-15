function nb = BoseFcn(x,T)
beta = 1/(KBoltzmann_ev*T);
nb = 1./(exp(x*beta)-1);
end