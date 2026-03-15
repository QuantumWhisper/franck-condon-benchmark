function M = generateMatrixW(rateW4DCache,N,vmode,alphaL,alphaR,lambda,Vsd,T,eta,Vg,tau)
% ,vmode,alpha,lambda,tau,T)
%{
generateMatrixW generates the matrix W as in 0 = WP, where W is a
coefficient matrix, P consists all steady state probabilities
In the first step, the steady-state occupation probabilities Pq
n are obtained from the rate equations 0 = dP^n_q/dt
N - vector space size, i.e q = 0,1,2,...,N-1
tau - time scale, equilibrated phonons tau→0;unequilibrated phonons tau→Inf
dependencies: rateW, rateW_lead
%}

M = zeros(2*N,2*N);

% for n = 0 case
base = 0;
for ii = 1:N
    m1Index = base + ii;
    for jj = 1:(2*N)
        m2Index = jj;
        modjj = mod(jj,N);
        if modjj == 0
            modjj = N;
        end
        P_q_index = modjj - 1;
        if jj <= N
            if ii == modjj
                n1 = 0;
                n2 = 0;
                q1 = P_q_index;
                q2 = P_q_index;
                M(m1Index,m2Index) = rateW_lead(rateW4DCache,N,n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,Vg)...
                    - sigmaW00(rateW4DCache,N,q1,vmode,alphaL,alphaR,lambda,Vsd,T,eta,Vg)...
                    - sigmaW01(rateW4DCache,N,q1,vmode,alphaL,alphaR,lambda,Vsd,T,eta,Vg)...
                    - 1/tau + peq_(P_q_index,vmode,T)/tau;
            else
                n1 = 0;
                n2 = 0;
                q1 = P_q_index;
                q2 = ii - 1;
                M(m1Index,m2Index) = rateW_lead(rateW4DCache,N,n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,Vg)...
                    + peq_(q2,vmode,T)/tau;
            end
        else
            n1 = 1;
            n2 = 0;
            q1 = P_q_index;
            q2 = ii - 1;
            M(m1Index,m2Index) = rateW_lead(rateW4DCache,N,n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,Vg);
        end
    end
end

% for n = 1 case
base = N;
for ii = 1:N
    m1Index = base + ii;
    for jj = 1:(2*N)
        modjj = mod(jj,N);
        if modjj == 0
            modjj = N;
        end
        if jj <= N
            n1 = 0;
            n2 = 1;
            q1 = modjj - 1;
            q2 = ii - 1;
            M(m1Index,jj) = rateW_lead(rateW4DCache,N,n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,Vg);
        else
            if modjj == ii
                n1 = 1;
                n2 = 1;
                q1 = modjj - 1;
                q2 = ii - 1;
                M(m1Index,jj) = rateW_lead(rateW4DCache,N,n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,Vg)...
                    -sigmaW10(rateW4DCache,N,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,Vg)...
                    -sigmaW11(rateW4DCache,N,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,Vg)...
                    -1/tau + peq_(q2,vmode,T)/tau;
            else
                n1 = 1;
                n2 = 1;
                q1 = modjj - 1;
                q2 = ii - 1;
                M(m1Index,jj) = rateW_lead(rateW4DCache,N,n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,Vg)...
                    + peq_(q2,vmode,T)/tau;
            end
        end       
    end
end

end

function tol = sigmaW11(rateW4DCache,N,q1,vmode,alphaL,alphaR,lambda,Vsd,T,eta,Vg)
% sum of q2 when n2 = 1
targetMatrix = zeros(1,N);
n1 = 1;
n2 = 1;
q2 = 0:(N-1);
targetMatrix(:) = rateW_lead(rateW4DCache,N,n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,Vg);
tol = sum(targetMatrix);
end

function tol = sigmaW10(rateW4DCache,N,q1,vmode,alphaL,alphaR,lambda,Vsd,T,eta,Vg)
% sum of q2 when n2 = 0
targetMatrix = zeros(1,N);
n1 = 1;
n2 = 0;
q2 = 0:(N-1);
targetMatrix(:) = rateW_lead(rateW4DCache,N,n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,Vg);
tol = sum(targetMatrix);
end

function p = peq_(q,vmode,T)

beta = 1/(KBoltzmann_ev*T);
p = exp(-q*vmode*beta).*(1-exp(-vmode*beta));

end

function tol = sigmaW00(rateW4DCache,N,q1,vmode,alphaL,alphaR,lambda,Vsd,T,eta,Vg)
% sum of q2 when n2 = 0
targetMatrix = zeros(1,N);
n1 = 0;
n2 = 0;
q2 = 0:(N-1);
targetMatrix(:) = rateW_lead(rateW4DCache,N,n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,Vg);
tol = sum(targetMatrix);
end

function tol = sigmaW01(rateW4DCache,N,q1,vmode,alphaL,alphaR,lambda,Vsd,T,eta,Vg)
% sum of q2 when n2 = 1
targetMatrix = zeros(1,N);
n1 = 0;
n2 = 1;
q2 = 0:(N-1);
targetMatrix(:) = rateW_lead(rateW4DCache,N,n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,Vg);
tol = sum(targetMatrix);
end
