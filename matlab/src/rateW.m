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
