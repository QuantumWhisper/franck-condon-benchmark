function I = regularizedI(E1,E2,epsilon1,epsilon2,T)
% [f,real_f,img_f] = digammaFcn(k,x)
beta = 1/(KBoltzmann_ev*T);
epsilon1 = epsilon1';
a1 = 1/2 + 1i*beta*(E2-epsilon1)/(2*pi);
a2 = 1/2 - 1i*beta*(E2-epsilon2)/(2*pi);
a3 = 1/2 + 1i*beta*(E1-epsilon1)/(2*pi);
a4 = 1/2 - 1i*beta*(E1-epsilon2)/(2*pi);
t1 = digammaFcn(0,a1);
t2 = digammaFcn(0,a2);
t3 = digammaFcn(0,a3);
t4 = digammaFcn(0,a4);
I = BoseFcn(E2-E1,T)./(epsilon1-epsilon2).*real(t1-t2-t3+t4);
end