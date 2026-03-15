function tol = sumMMr11(q1,q2,lambda,muL,muR,vmode,epsilond,T)
tempN = round(lambda^2.2*3);
epsilon = 1e-14;
temp_tol = m_sumMMr11(tempN,q1,q2,lambda,muL,muR,vmode,epsilond,T);
tempN = tempN + 5;
temp_tol2 = m_sumMMr11(tempN,q1,q2,lambda,muL,muR,vmode,epsilond,T);
relative_diff = abs(log10(temp_tol2./temp_tol));
while any(relative_diff > epsilon)
    locs = relative_diff > epsilon;
    tempN = tempN + min(max(round(tempN*0.5),10),20);
    temp_tol = temp_tol2;
    temp_tol2_temp = m_sumMMr11(tempN,q1,q2(locs),lambda,muL,muR,vmode,epsilond,T);
    temp_tol2(locs) = temp_tol2_temp;
    relative_diff = abs(log10(temp_tol2./temp_tol));
end
%fprintf('sumMMr11(q1 = %g, nq2 = %g): Final tempN = %g with relative difference (max: %.4g) %.4g +- %.4g\n',q1, numel(q2),tempN,max(relative_diff),mean(relative_diff,'omitmissing'),std(relative_diff,'omitmissing'))
tol = temp_tol2;
end

function result = temp_sumMMr11(varargin)

temp = memoize(@m_sumMMr11);
temp.CacheSize = 1e+8;
result = temp(varargin{:});

end

