function [temp] = test_sumMMr
N = 10;
vmode = 73e-3;
alphaL = 0.02;
alphaR = 0.02;
lambda = 5;
T = 4.2;
eta = 1/2;
lead = 1; % -1
Vg = 0.1;
Vsd = 0.1;
gammaL = alphaL*vmode;
gammaR = alphaR*vmode;
epsilond = 0 + Vg;
muL = eta*Vsd;
muR = -(1-eta)*Vsd;
q1 = 4;
q2 = 0:9;
xdata2plot = 1:1:numel(q2);
ydata2plot = zeros(size(xdata2plot));
nx = numel(xdata2plot);
temp = cell([1 nx]);
time2plot = zeros(size(xdata2plot));
mFcn_sumMMMMrs = memoize(@sumMMMMrs);
mFcn_sumMMMMrs.CacheSize = 1e+3;
temp_sumMMMMrs = memoize(@m_sumMMMMrs);
temp_sumMMMMrs.CacheSize = 1e+6;
temp_m_sumMMr = memoize(@m_sumMMr);
temp_m_sumMMr.CacheSize = 1e+6;
temp_m_sumMMr11 = memoize(@m_sumMMr11);
temp_m_sumMMr11.CacheSize = 1e+6;
temp_sumMMMMrs11 = memoize(@m_sumMMMMrs11);
temp_sumMMMMrs11.CacheSize = 1e+6;
temp_sumMMr = memoize(@sumMMr);
temp_sumMMr.CacheSize = 1e+6;
temp_sumMMr11 = memoize(@sumMMr11);
temp_sumMMr11.CacheSize = 1e+6;
tic

tol = temp_sumMMr11(q1,q2,lambda,muL,muR,vmode,epsilond,T);
ydata2plot = tol;

runtime = toc;
fprintf('Total runtime: %.3g s\n',runtime)

fig = CMDA_createTagedFigure('test_sumMMr');
cax1 = subplottight(1,1,1);
plot(cax1,xdata2plot,ydata2plot)

end