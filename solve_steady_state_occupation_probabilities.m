function P = solve_steady_state_occupation_probabilities(M)

%{
solves the linear system C*x = d in the least-squares sense, subject to A*x ≤ b.
x = lsqlin(C,d,A,b,Aeq,beq,lb,ub) adds linear equality constraints Aeq*x = beq 
and bounds lb ≤ x ≤ ub. If you do not need certain constraints such as Aeq and beq, 
set them to []. If x(i) is unbounded below, set lb(i) = -Inf, 
and if x(i) is unbounded above, set ub(i) = Inf.
P - steady state probabilities
P^0_0
P^0_1
...
P^0_i
...
P^0_(N-1)
P^1_0
P^1_1
...
P^1_i
...
P^1_(N-1)
size of 2N
%}
[~,N] = size(M);
d = zeros(1,N);
C = M;
A = [];
b = [];
beq = 1;
Aeq = ones(1,N);
lb = zeros(N,1);
ub = ones(N,1);
x0 = [];
%options = struct;
options = optimoptions('lsqlin','Algorithm','interior-point','MaxIterations',2e+3);
P = lsqlin(C,d,A,b,Aeq,beq,lb,ub,[],options);

end