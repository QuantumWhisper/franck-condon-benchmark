function targetMatrixQ2 = sumMMMMrs(q1,q2,lambda,muL,muR,vmode,epsilond,T)
nq2 = numel(q2);
Q2 = q2;
targetMatrixQ2 = zeros(1,nq2);
temp_sumMMMMrs = memoize(@m_sumMMMMrs);
temp_sumMMMMrs.CacheSize = 1e+6;
parfor qq2 = 1:nq2    
    tempN = round(lambda^2*4);
    q2 = Q2(qq2);
    epsilon = 1e-14;
    temp_tol = m_sumMMMMrs(tempN,q1,q2,lambda,muL,muR,vmode,epsilond,T);
    tempN = tempN + 5;
    temp_tol2 = m_sumMMMMrs(tempN,q1,q2,lambda,muL,muR,vmode,epsilond,T);
    relative_diff = abs(log(temp_tol2/temp_tol));
    while relative_diff > epsilon
        temp_tol = temp_tol2;
        tempN = tempN + min(max(round(tempN*0.5),20),40);
        temp_tol2 = m_sumMMMMrs(tempN,q1,q2,lambda,muL,muR,vmode,epsilond,T);
        relative_diff = abs(log(temp_tol2/temp_tol));
    end
    %fprintf('sumMMMMrs(q1 = %g q2 = %g): tempN = %g with relative difference %.4g\n',q1, q2, tempN, relative_diff)
    tol = temp_tol2;
    targetMatrixQ2(qq2) = tol;
end

end


