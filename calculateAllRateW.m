function [rateW4D,fpath] = calculateAllRateW(N,vmode,alphaL,alphaR,lambda,Vsd,T,eta,lead,Vg)

% check it if has been calculated before
fname = 'rateWcache.mat';
folder = pwd;
fpath = fullfile(folder,fname);
if exist(fpath,"file") == 2
    temp = load(fpath);
    temp_m_rateW = temp.temp_m_rateW;
    temp_m_rateW_stats = temp_m_rateW.stats();
    if temp_m_rateW_stats.CacheOccupancyPercent > 0.80
        warning('Cache in rateWcache.mat is 80% full.')
        temp_m_rateW.CacheSize = 2 * temp_m_rateW.CacheSize;
    end
else

end

% for each dimention q1,q2,n1,n2
rateW4D = zeros([N N 2 2]);
for n1 = 0:1
    for n2 = 0:1
        targetMatrix = zeros(N,N);
        for ii = 1:N
            q1 = ii - 1;
            q2 = 0:(N-1);
            targetMatrix(ii,:) = rateW([],n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,lead,Vg);
        end
        rateW4D(:,:,n1+1,n2+1) = targetMatrix;
    end
end

assignin("base","rateW4DCache",rateW4D)
temp_m_rateW = memoize(@m_rateW);
save(fpath,'temp_m_rateW')
fprintf('\nCache successfully saved at:%s\n',fpath)
end