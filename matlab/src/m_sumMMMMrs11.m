function temp_tol = m_sumMMMMrs11(N,q1,q2,lambda,muL,muR,vmode,epsilond,T)
rr = 0:N-1;
ss = 0:(N-1);
%MMMM_quad = FCMatrix(q2,rr,lambda).*conj(FCMatrix(q1,rr,lambda)).*conj(FCMatrix(q2,ss,lambda)).*FCMatrix(q1,ss,lambda);
MMMM_quad2D = FCMatrix(q2,rr,lambda)'.*conj(FCMatrix(q1,rr,lambda))'.*conj(FCMatrix(q2,ss,lambda)).*FCMatrix(q1,ss,lambda);
Irs = regularizedI(muL,muR-(q1-q2)*vmode,epsilond+(q2-rr)*vmode,epsilond+(q2-ss)*vmode,T);
targetMatrix = MMMM_quad2D.*Irs;
targetMatrix(1:1+size(targetMatrix,1):end) = 0; % very powerful, all diag elements equal to 0
temp_tol = sum(targetMatrix,"all");
end