function w = rateW(varargin)
%{
rateW calculates the rate for changing the charge number of the molecule 
from n1 to n2 by tunneling across lead, and simultaneously
changing the number of excited phonons from q1 to q2
N - 0, 1, ... N-1 phonon states
vmode - \hbar\omega phonon mode energy in eV, typical value 0.073 eV
alphaL/alphaR - ratio of tunnel width of left or right lead and vmode,
typical value 0.02 as in https://journals.aps.org/prb/pdf/10.1103/PhysRevB.74.205438
lambda - e-p coupling strength
Vsd - source drain bias, V
T - temperature, K
eta - symmetric junction eta = 1/2 corresponds to ul = -ur = Vsd/2
lead
Vg - gate
see: generateMatrixW, rateW_lead
%}
% check it if has been calculated before
%{
fname = 'rateW4DCache.mat';
folder = pwd;
fpath = fullfile(folder,fname);
rateW4DON = false;
if exist(fpath,"file") == 2
    temp = load(fpath);
    rateW4DCache = temp.rateW4DCache;
else
    rateW4DCache = struct;
end
value = [N,vmode,alphaL,alphaR,lambda,Vsd,T,eta,lead,Vg];
[rateW4DCache,settingValue,flag] = get_set_values(rateW4DCache,'setting',value);
if ~flag

else
    rateW4D = rateW4DCache.rateW4Ds{settingValue};
    rateW4DON = true;
end
if rateW4DON
    w = rateWFast(rateW4D,n1,n2,q1,q2);
    return
end
%}
%{
if isempty(rateW4DCache)
else
    value = [N,vmode,alphaL,alphaR,lambda,Vsd,T,eta,lead,Vg];
    [rateW4DCache,settingValue,~] = get_set_values(rateW4DCache,'setting',value);
    rateW4D = rateW4DCache.rateW4Ds{settingValue};
    w = rateWFast(rateW4D,n1,n2,q1,q2);
    return
end
%}
%fprintf('seeing this...\n')
%T = 4.2; % K
%N = max([q1,q2]);
%N = min(max(50,N*4),N*10); % used for sum convergence 

rateW4DCache = varargin{1};
if isempty(rateW4DCache)
    temp = memoize(@m_rateW);
    if temp.CacheSize <= 1e+10
        temp.CacheSize = 1e+10;
    end 
    w = temp(varargin{2:end});
else
    w = rateW4DCache(varargin{2:end});
end

end

