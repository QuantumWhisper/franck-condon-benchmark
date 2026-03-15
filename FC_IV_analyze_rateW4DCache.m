function [vmode,alphaL,alphaR,lambda,Vsd,T,eta,lead,Vg] = FC_IV_analyze_rateW4DCache(varargin)

fpath = fullfile(pwd,'rateWcache.mat');
if numel(varargin) > 0
    if exist(varargin{1},"file") == 2
        fpath = varargin{1};
    end
end

temp = load(fpath);
temp_m_rateW = temp.temp_m_rateW;
temp_m_rateW_stats = temp_m_rateW.stats();
temp_parametersMatrix = cellfun(@(x) cell2mat(x),temp_m_rateW_stats.Cache.Inputs','UniformOutput',false);
parametersMatrix = cell2mat(temp_parametersMatrix);
% ref: value = (N,n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,lead,Vg); in
% m_rateW.m

vmode = unique(parametersMatrix(:,end-8));
alphaL = unique(parametersMatrix(:,end-7));
alphaR = unique(parametersMatrix(:,end-6));
lambda = unique(parametersMatrix(:,end-5));
Vsd = unique(parametersMatrix(:,end-4));
T = unique(parametersMatrix(:,end-3));
eta = unique(parametersMatrix(:,end-2));
lead = unique(parametersMatrix(:,end-1));
Vg = unique(parametersMatrix(:,end));

end