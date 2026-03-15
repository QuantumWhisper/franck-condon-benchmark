function temp_tol = m_sumMMMMrs(N,q1,q2,lambda,muL,muR,vmode,epsilond,T)

rr = 0:N-1;
ss = 0:(N-1);
targetMatrix = MMMM_quad_Irs(q1,q2,lambda,muL,muR,vmode,epsilond,T,rr,ss);
targetMatrix(1:1+size(targetMatrix,1):end) = 0; % very powerful, all diag elements equal to 0
temp_tol = sum(targetMatrix,"all");


end

function result = MMMM_quad_Irs(varargin)

temp = memoize(@m_MMMM_quad_Irs);
temp.CacheSize = 1e+6;
result = temp(varargin{:});

end