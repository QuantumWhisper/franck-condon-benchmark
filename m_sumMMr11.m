function temp_tol = m_sumMMr11(N,q1,q2,lambda,muL,muR,vmode,epsilond,T)
rr = 0:N-1;
%MM_square = (FCMatrix(q2,rr,lambda).*FCMatrix(q1,rr,lambda)).^2;
MM_square2D = abs(FCMatrix(q2,rr,lambda).*FCMatrix(q1,rr,lambda)).^2';
%Jr = regularizedJ11(muL,muR-(q1-q2)*vmode,epsilond+q2_rr*vmode,T);
Jr2D = regularizedJ(muL,muR-(q1-q2)*vmode,epsilond+(q2'-rr)*vmode,T);
targetMatrix = MM_square2D.*Jr2D;
temp_tol = sum(targetMatrix,1);
end