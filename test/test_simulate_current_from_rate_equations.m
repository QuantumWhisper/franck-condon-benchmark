function test_simulate_current_from_rate_equations
%{
% test rateW
n1 = 0;
n2 = 0;
q1 = 0;
q2 = 1;
vmode = 73e-3;
alphaL = 0.02;
alphaR = 0.02;
lambda = 5;
Vsd = 0;
T = 4.2; 
eta = 1/2;
lead = 1; % -1
Vg = 0;
N = 30;
w = rateW(N,n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,lead,Vg)
%}
%{
% test rateW
n1 = 0;
n2 = 1;
q1 = 4;
q2 = 3;
vmode = 73e-3;
alphaL = 0.02;
alphaR = 0.02;
lambda = 5;
Vsd = 0;
T = 4.2; 
eta = 1/2;
lead = 1; % -1
Vg = 0;
N = 5;
w = rateW_lead(N,n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,lead,Vg)
%}
%{
% test M = generateMatrixW(N,vmode,alphaL,alphaR,lambda,Vsd,T,eta,lead,Vg,tau)
N = 2;
vmode = 73e-3;
alphaL = 0.02;
alphaR = 0.02;
lambda = 5;
Vsd = 0.1;
T = 4.2; 
eta = 1/2;
lead = 1; % -1
Vg = 0;
tau = 0.001;
M = generateMatrixW(N,vmode,alphaL,alphaR,lambda,Vsd,T,eta,lead,Vg,tau)
P = solve_steady_state_occupation_probabilities(M)
%}
% test [I_tol,I_seq,I_cot] = current_from_rate_equations(Vsd,N,vmode,alphaL,alphaR,lambda,T,eta,lead,Vg,tau)
%{
Vsd = 0;
N = 2;
vmode = 73e-3;
alphaL = 0.02;
alphaR = 0.02;
lambda = 5;
%Vsd = 0;
T = 4.2; 
eta = 1/2;
lead = 1; % -1
Vg = 0;
tau = Inf;
tic
[I_tol,I_seq,I_cot] = current_from_rate_equations(Vsd,N,vmode,alphaL,alphaR,lambda,T,eta,lead,Vg,tau)
toc
%}
end