function temp_tol = m_sumMMr(N,q1,q2,lambda,muL,muR,vmode,epsilond,T)
%nq2 = numel(q2);
%targetMatrix = zeros(N,nq2);
rr = 0:N-1;
%MM_square = abs(FCMatrix(q2,rr,lambda).*FCMatrix(q1,rr,lambda)).^2;
MM_square2D = abs(FCMatrix(q2,rr,lambda)'.*FCMatrix(q1,rr,lambda)').^2;
Jr = regularizedJ(muL,muR-(q1-q2)*vmode,epsilond-(q1-rr)*vmode,T);
targetMatrix = MM_square2D.*Jr;
temp_tol = sum(targetMatrix,1);
end