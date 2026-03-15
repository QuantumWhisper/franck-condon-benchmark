function [I_tol,I_seq,I_cot] = current_from_rate_equations(Vsd,N,vmode,alphaL,alphaR,lambda,T,eta,lead,Vg,tau,rateW4DCache)

M = generateMatrixW(rateW4DCache,N,vmode,alphaL,alphaR,lambda,Vsd,T,eta,lead,Vg,tau);
P = solve_steady_state_occupation_probabilities(M);

leadR = -1;
leadL = 1;
temp = 0;
% for n = 0 case
targetMatrix = zeros(N,N);
for ii = 1:N
    n1 = 0;
    n2 = 1;
    q1 = ii - 1;
    q2 = 0:(N-1);
    w01R = rateW(rateW4DCache,n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,leadR,Vg);
    w01L = rateW(rateW4DCache,n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,leadL,Vg);
    diffW01 = w01R - w01L;
    targetMatrix(ii,:) = psub_(P,0,q1).*diffW01;
end
I_seq0 = sum(targetMatrix,"all");
%{
for q1 = 0:(N-1)
    for q2 = 0:(N-1)
        n1 = 0;
        n2 = 1;
        w01R = rateW(rateW4DCache,N,n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,leadR,Vg);
        w01L = rateW(rateW4DCache,N,n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,leadL,Vg);
        diffW01 = w01R - w01L;
        temp = temp + psub_(P,0,q1)*diffW01;
    end
end
I_seq0 = temp;
%}
% for n = 1 case
temp = 0;
targetMatrix = zeros(N,N);
for ii = 1:N
    n1 = 1;
    n2 = 0;
    q1 = ii - 1;
    q2 = 0:(N-1);
    w10R = rateW(rateW4DCache,n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,leadR,Vg);
    w10L = rateW(rateW4DCache,n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,leadL,Vg);
    diffW01 = w10L - w10R;
    targetMatrix(ii,:) = psub_(P,1,q1).*diffW01;
end
I_seq1 = sum(targetMatrix,"all");
%{
for q1 = 0:(N-1)
    for q2 = 0:(N-1)
        n1 = 1;
        n2 = 0;
        w10R = rateW(rateW4DCache,N,n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,leadR,Vg);
        w10L = rateW(rateW4DCache,N,n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,leadL,Vg);
        diffW01 = w10L - w10R;
        temp = temp + psub_(P,1,q1)*diffW01;
    end
end
I_seq1 = temp;
%}
I_seq = I_seq0 + I_seq1;
% for cotunneling case
temp = 0;
for nn = 0:1
    %{
    targetMatrix = zeros(N,N);
    parfor ii = 1:N
        q1 = ii - 1;
        q2 = 0:(N-1);
        p_nn_q1 = psub_(P,nn,q1);
        n1 = nn;
        n2 = nn;        
        wRL = rateW(rateW4DCache,N,n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,leadR,Vg);
        wLR = rateW(rateW4DCache,N,n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,leadL,Vg);
        diffW_nn = wRL - wLR;
        targetMatrix(ii,:) = p_nn_q1.*diffW_nn;
    end
    temp = temp + sum(targetMatrix,"all");
    %}
    for q1 = 0:(N-1)
        q2 = 0:(N-1);
        p_nn_q1 = psub_(P,nn,q1);
        n1 = nn;
        n2 = nn;
        wRL = rateW(rateW4DCache,n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,leadR,Vg);
        wLR = rateW(rateW4DCache,n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,leadL,Vg);
        diffW_nn = wRL - wLR;
        temp = temp + sum(p_nn_q1*diffW_nn);
    end
    %}
end
I_cot = temp;
I_tol = I_seq + I_cot;

end
function ps = psub_(P,n,q)

N = numel(P)/2;
if n == 0
    base = 0;
else
    base = N;
end
qIndex = q + 1 + base;
ps = P(qIndex);


end