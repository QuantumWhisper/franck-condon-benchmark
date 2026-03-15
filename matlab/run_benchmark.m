function results = run_benchmark()
% RUN_BENCHMARK Run Franck-Condon I-V simulation benchmark
%   Loads parameters from benchmark/spec/default_params.json,
%   runs the simulation with timing, and saves results as JSON.

% Setup paths
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, 'src'));
addpath(fullfile(thisDir, 'utils'));

% Load benchmark parameters
specPath = fullfile(thisDir, '..', 'benchmark', 'spec', 'default_params.json');
spec = jsondecode(fileread(specPath));
params = spec.parameters;
sweep = spec.bias_sweep;

% Extract parameters
N = params.N;
vmode = params.vmode;
alphaL = params.alphaL;
alphaR = params.alphaR;
lambda = params.lambda;
T = params.T;
eta = params.eta;
Vg = params.Vg;
if ischar(params.tau) || isstring(params.tau)
    tau = Inf;
else
    tau = params.tau;
end

% Build bias sweep
Vsd = sweep.Vsd_start:sweep.Vsd_step:sweep.Vsd_end;
nVsd = numel(Vsd);

% Pre-allocate output arrays
I_tol = zeros(size(Vsd));
I_seq = zeros(size(Vsd));
I_cot = zeros(size(Vsd));

% --- Begin timed simulation ---
fprintf('=== Franck-Condon Benchmark (MATLAB) ===\n');
fprintf('N=%d, lambda=%.1f, T=%.1f K, Vsd=[%.3f:%.3f:%.3f] V\n', ...
    N, lambda, T, sweep.Vsd_start, sweep.Vsd_step, sweep.Vsd_end);
fprintf('Total bias points: %d\n\n', nVsd);

tic;
for vv = 1:nVsd
    v = Vsd(vv);
    fprintf('Vsd = %.4f V  (%d/%d)\n', v, vv, nVsd);

    % Pre-compute rate matrices for both leads
    calculateAllRateW(N, vmode, alphaL, alphaR, lambda, v, T, eta, 1, Vg);
    calculateAllRateW(N, vmode, alphaL, alphaR, lambda, v, T, eta, -1, Vg);
    rateW4DCache = [];

    % Compute current
    [I_tol(vv), I_seq(vv), I_cot(vv)] = current_from_rate_equations( ...
        v, N, vmode, alphaL, alphaR, lambda, T, eta, Vg, tau, rateW4DCache);

    % Sign convention: current flows opposite to bias
    I_tol(vv) = -sign(v) * I_tol(vv);
    I_seq(vv) = -sign(v) * I_seq(vv);
    I_cot(vv) = -sign(v) * I_cot(vv);
end
wall_time = toc;
% --- End timed simulation ---

% Convert to SI (Amperes)
e = ee_ElementaryCharge();
I_tol = e * I_tol;
I_seq = e * I_seq;
I_cot = e * I_cot;

% Build results struct
results = struct();
results.language = 'MATLAB';
results.version = version();
results.wall_time_seconds = wall_time;
results.parameters = params;
results.bias_sweep = sweep;
results.Vsd = Vsd;
results.I_tol = I_tol;
results.I_seq = I_seq;
results.I_cot = I_cot;
results.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');

% Save results as JSON
resultsDir = fullfile(thisDir, '..', 'benchmark', 'results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end
resultsPath = fullfile(resultsDir, 'matlab_results.json');
fid = fopen(resultsPath, 'w');
fprintf(fid, '%s', jsonencode(results, 'PrettyPrint', true));
fclose(fid);

% Print summary
fprintf('\n=== Benchmark Complete ===\n');
fprintf('Wall time: %.2f seconds\n', wall_time);
fprintf('Results saved to: %s\n', resultsPath);
fprintf('Max |I_tol|: %.4e A\n', max(abs(I_tol)));
fprintf('Max |I_seq|: %.4e A\n', max(abs(I_seq)));
fprintf('Max |I_cot|: %.4e A\n', max(abs(I_cot)));

end
