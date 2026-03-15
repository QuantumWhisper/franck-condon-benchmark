function rateW4D = test_calculateAllRateW
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

tic

rateW4D = zeros([N N 2 2]);
for n1 = 0:1
    for n2 = 0:1
        targetMatrix = zeros(N,N);
        for ii = 1:N
            q1 = ii - 1;
            q2 = 0:(N-1);
            targetMatrix(ii,:) = rateW([],n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,lead,Vg);
        end
        rateW4D(:,:,n1+1,n2+1) = targetMatrix;
    end
end

runtime = toc;
fprintf('Total runtime: %.3g s\n',runtime)

fig = CMDA_createTagedFigure('test_calculateAllRateW');
cax1 = subplottight(1,1,1);

q1 = 1:1:N;
q2 = 1:1:N;
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


end