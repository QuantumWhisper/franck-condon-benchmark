function fig = simulate_current_from_rate_equations
% simulate IV characteristics for Franck-Condon regime
% based on quantum master equations
% last update date: 2023-May-21

tic
[N,vmode,alphaL,alphaR,lambda,Vsd,T,eta,lead,Vg,tau] = setParameters;
Vsd = 0:0.003:0.6;
N = 6;
lambda = 5;
fig = CMDA_createTagedFigure(': simulation of IV from rate equations');
cax1 = subplottight(1,1,1);
nVsd = numel(Vsd);

I_tol = zeros(size(Vsd));
I_seq = zeros(size(Vsd));
I_cot = zeros(size(Vsd));


for vv = 1:nVsd    
    v = Vsd(vv);
    fprintf('Working bias voltage = %g V\n',v)
    %[rateW4D,fpath] = 
    calculateAllRateW(N,vmode,alphaL,alphaR,lambda,v,T,eta,1,Vg);
    calculateAllRateW(N,vmode,alphaL,alphaR,lambda,v,T,eta,-1,Vg);
    rateW4DCache = [];
    [I_tol(vv),I_seq(vv),I_cot(vv)] = current_from_rate_equations(v,N,vmode,alphaL,alphaR,lambda,T,eta,lead,Vg,tau,rateW4DCache);
    I_tol(vv) = -sign(v)*I_tol(vv);
    I_seq(vv) = -sign(v)*I_seq(vv);
    I_cot(vv) = -sign(v)*I_cot(vv);
    fprintf('running time %g s\n',toc)
end
I_tol = ee_ElementaryCharge*I_tol;
I_seq = I_seq*ee_ElementaryCharge;
I_cot = I_cot*ee_ElementaryCharge;
assignin("base","I_tol",I_tol)
assignin("base","I_seq",I_seq)
assignin("base","I_cot",I_cot)
plot(cax1,Vsd,I_tol,'o')
hold(cax1,"on")
plot(cax1,Vsd,I_seq)
hold(cax1,"on")
plot(cax1,Vsd,I_cot)
lgr = {'I_{tol}','I_{seq}','I_{cot}'};
legend(cax1,lgr)
titlestr = sprintf('N = %g $\\hbar\\omega$ = %g V $\\lambda$ = %g T = %g K $\\tau$ = %g $\\frac{\\Gamma_{L}}{\\hbar\\omega}$ = %g $\\frac{\\Gamma_{R}}{\\hbar\\omega}$ = %g $\\eta$ = %g $V_g$ = %g',...
    N,vmode,lambda,T,tau,alphaL,alphaR,eta,Vg);
title(cax1,titlestr,'interpreter','latex');
xlabel(cax1,'Voltage (V)')
ylabel(cax1,'Current (A)')
fprintf('Finish running simulate_current_from_rate_equations in %.2f s\n',toc)

end
function [N,vmode,alphaL,alphaR,lambda,Vsd,T,eta,lead,Vg,tau] = setParameters
tau = Inf;%;Inf;%1e-9;

try
    [vmode,alphaL,alphaR,lambda,Vsd,T,eta,lead,Vg] = FC_IV_analyze_rateW4DCache;
    % [N,vmode,alphaL,alphaR,lambda,Vsd,T,eta,lead,Vg]
    Vg = 0;% Vg(1);
    N = 8;
    %error('0')
catch ME
    fprintf('error while loading rateW4DCache.mat in pwd.\nDefault will be used.\n')
    %rethrow(ME)
    %return
    N = 8;
    vmode = 73e-3;
    alphaL = 0.02;
    alphaR = 0.02;
    lambda = 1.5;
    T = 4.2;
    eta = 1/2;
    lead = 1; % -1
    Vg = 0;
    Vsd = 0:0.01:0.6;%:0.02:0.3;
end
end