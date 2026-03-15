function [rateW4D,temp] = test_m_sumMMr
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
q1 = 1;
q2 = 0:9;
xdata2plot = 1:1:100;
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
tic

for ii = 1:nx
    %tol = mFcn_sumMMMMrs(ii,q1,q2,lambda,muL,muR,vmode,epsilond,T);
    %q2 =3;
    %tol = temp_sumMMMMrs(ii,q1,q2,lambda,muL,muR,vmode,epsilond,T);
    tol = temp_m_sumMMr(ii,q1,q2,lambda,muL,muR,vmode,epsilond,T);
    %tol = sumMMr11(N,q1,q2,lambda,muL,muR,vmode,epsilond,T);
    temp{ii} = tol;
    if ii == 1
        relative_diff = 0;
    else
        temp_tol2 = tol(1);
        temp_tol = temp{ii-1}(1);
        relative_diff = abs(log(temp_tol2/temp_tol)/log(temp_tol));
    end
    ydata2plot(ii) = relative_diff;
    runtime = toc;
    time2plot(ii) = runtime;
    fprintf('ii = %d runtime: %.3g s\n',ii,runtime)
end

rateW4D = zeros([N N 2 2]);
%{
for n1 = 0:1
    for n2 = 0:1
        targetMatrix = zeros(N,N);
        for ii = 1:N
            q1 = ii - 1;
            %q2 = 0:(N-1);
            %targetMatrix(ii,:) = rateW([],N,n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,lead,Vg);
            for qq = 0:(N-1)
                w = rateW([],N,n1,n2,q1,qq,vmode,alphaL,alphaR,lambda,Vsd,T,eta,lead,Vg);
            end
            runtime = toc;
            fprintf('n1 = %g n2 = %g q1 = %g running time: %.3g s\n',n1,n2,q1,runtime)
        end
        rateW4D(:,:,n1+1,n2+1) = targetMatrix;
    end
end
%}
runtime = toc;
fprintf('Total runtime: %.3g s\n',runtime)

fig = CMDA_createTagedFigure('temp');
cax1 = subplottight(1,1,1);
plot(cax1,xdata2plot,ydata2plot)
%{
q1 = 1:1:10;
q2 = 1:1:10;
[Q2,Q1] = meshgrid(q2,q1);
surfaceM = rateW4D(:,:,2,2);
s = surface(cax1,Q1,Q2,surfaceM);
%s.EdgeColor = 'none';
xlabel(cax1,'q_1')
ylabel(cax1,'q_2')
zlabel(cax1,'w_{q_1q_2}')
title(cax1,sprintf('\\lambda = %g',lambda))
hcb=colorbar;
hcb.Title.String = 'w_{q_1q_2}';
%}


end