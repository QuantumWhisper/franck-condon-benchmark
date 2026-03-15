function result = m_MMMM_quad_Irs(q1,q2,lambda,muL,muR,vmode,epsilond,T,rr,ss)

%MMMM_quad = FCMatrix(q2,rr,lambda).*conj(FCMatrix(q1,rr,lambda)).*conj(FCMatrix(q2,ss,lambda)).*FCMatrix(q1,ss,lambda);
MMMM_quad2D = FCMatrix(q2,rr,lambda)'.*conj(FCMatrix(q1,rr,lambda))'.*conj(FCMatrix(q2,ss,lambda)).*FCMatrix(q1,ss,lambda);
Irs = regularizedI(muL,muR-(q1-q2)*vmode,epsilond-(q1-rr)*vmode,epsilond-(q1-ss)*vmode,T);
result = MMMM_quad2D.*Irs;

end