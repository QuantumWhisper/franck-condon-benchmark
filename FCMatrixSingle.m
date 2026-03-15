function M = FCMatrixSingle(q1,q2,lambda)

q = min(q1,q2);
Q = max(q1,q2);
% laguerreL(n,x) and laguerreL(n,0,x) are equivalent.

l = laguerreL(q,Q-q,lambda.^2);

M = (sign(q2-q1)).^(q1-q2).*lambda.^(Q-q).*exp(-lambda.^2/2).*sqrt(factorial(q)./factorial(Q)).*l;

end