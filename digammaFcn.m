function [f,real_f,img_f] = digammaFcn(k,x)

symx = sym(x);
f = double(psi(k,symx));
real_f = real(f);
img_f = imag(f);

end