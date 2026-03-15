function tol = test_sumMMMMrs
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
q1 = 6;
q2 = 0:9;
xdata2plot = 1:1:150;
ydata2plot = zeros(size(xdata2plot));
nx = numel(xdata2plot);
temp = cell([1 nx]);
time2plot = zeros(size(xdata2plot));
mFcn_sumMMMMrs = memoize(@sumMMMMrs);
mFcn_sumMMMMrs.CacheSize = 1e+6;
mFcn_sumMMMMrs11 = memoize(@sumMMMMrs11);
mFcn_sumMMMMrs11.CacheSize = 1e+6;
temp_sumMMMMrs = memoize(@m_sumMMMMrs);
temp_sumMMMMrs.CacheSize = 1e+6;
temp_m_sumMMr = memoize(@m_sumMMr);
temp_m_sumMMr.CacheSize = 1e+6;
temp_m_sumMMr11 = memoize(@m_sumMMr11);
temp_m_sumMMr11.CacheSize = 1e+6;
temp_sumMMMMrs11 = memoize(@m_sumMMMMrs11);
temp_sumMMMMrs11.CacheSize = 1e+6;
tic
tol = mFcn_sumMMMMrs(q1,q2,lambda,muL,muR,vmode,epsilond,T);
runtime = toc;
fprintf('Total runtime: %.3g s\n',runtime)
fig = CMDA_createTagedFigure('test_m_sumMMMMrs');
cax1 = subplottight(1,1,1);
xdata2plot = 1:numel(tol);
ydata2plot = tol;
plot(cax1,xdata2plot,ydata2plot)

end