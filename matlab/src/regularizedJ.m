function J = regularizedJ(E1,E2,epsilon,T)
beta = 1/(KBoltzmann_ev*T);
epsilon = epsilon';
a1 = 1/2 + 1i*beta*(E2-epsilon)/(2*pi);
a2 = 1/2 + 1i*beta*(E1-epsilon)/(2*pi);
t1 = digammaFcn(1,a1);
t2 = digammaFcn(1,a2);
J = beta/(2*pi)*BoseFcn(E2-E1,T).*imag(t1-t2);
end