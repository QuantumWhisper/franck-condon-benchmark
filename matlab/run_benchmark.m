function results = run_benchmark(spec_name)
% RUN_BENCHMARK Run Franck-Condon I-V simulation benchmark
%   results = run_benchmark()          uses 'default' spec (N=15, publication quality)
%   results = run_benchmark('quick')   uses 'quick' spec (N=6, fast validation)
%   results = run_benchmark('default') uses 'default' spec explicitly
%
%   Outputs:
%     benchmark/results/matlab_<spec>_results.json   — full results
%     benchmark/results/matlab_<spec>_IV.csv          — I-V data for reuse
%     benchmark/results/matlab_<spec>_IV.pdf           — publication-quality plot (vector)
%     benchmark/results/matlab_<spec>_IV.png           — publication-quality plot (300 dpi)

if nargin < 1
    spec_name = 'default';
end

% Validate spec name
valid_specs = {'default', 'quick'};
if ~ismember(spec_name, valid_specs)
    error('Unknown spec "%s". Use "default" or "quick".', spec_name);
end

% Setup paths
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, 'src'));
addpath(fullfile(thisDir, 'utils'));

% Load benchmark parameters
specFile = sprintf('%s_params.json', spec_name);
specPath = fullfile(thisDir, '..', 'benchmark', 'spec', specFile);
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
fprintf('=== Franck-Condon Benchmark (MATLAB) [%s] ===\n', spec_name);
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
results.spec = spec_name;
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

% --- Save results ---
resultsDir = fullfile(thisDir, '..', 'benchmark', 'results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end
prefix = sprintf('matlab_%s', spec_name);

% 1. JSON
jsonPath = fullfile(resultsDir, [prefix '_results.json']);
fid = fopen(jsonPath, 'w');
fprintf(fid, '%s', jsonencode(results, 'PrettyPrint', true));
fclose(fid);

% 2. CSV (for easy reuse in other tools)
csvPath = fullfile(resultsDir, [prefix '_IV.csv']);
csvData = [Vsd(:), I_tol(:), I_seq(:), I_cot(:)];
fid = fopen(csvPath, 'w');
fprintf(fid, 'Vsd_V,I_tol_A,I_seq_A,I_cot_A\n');
fprintf(fid, '%.6e,%.6e,%.6e,%.6e\n', csvData');
fclose(fid);

% 3. Publication-quality I-V plot (PDF + PNG)
pdfPath = fullfile(resultsDir, [prefix '_IV.pdf']);
pngPath = fullfile(resultsDir, [prefix '_IV.png']);
plot_IV(Vsd, I_tol, I_seq, I_cot, N, lambda, T, vmode, alphaL, alphaR, eta, Vg, ...
    wall_time, pdfPath, pngPath);

% Print summary
fprintf('\n=== Benchmark Complete [%s] ===\n', spec_name);
fprintf('Wall time:  %.2f seconds\n', wall_time);
fprintf('Max |I_tol|: %.4e A\n', max(abs(I_tol)));
fprintf('Max |I_seq|: %.4e A\n', max(abs(I_seq)));
fprintf('Max |I_cot|: %.4e A\n', max(abs(I_cot)));
fprintf('\nSaved:\n');
fprintf('  JSON: %s\n', jsonPath);
fprintf('  CSV:  %s\n', csvPath);
fprintf('  PDF:  %s\n', pdfPath);
fprintf('  PNG:  %s\n', pngPath);

end

function plot_IV(Vsd, I_tol, I_seq, I_cot, N, lambda, T, vmode, alphaL, alphaR, eta, Vg, ...
    wall_time, pdfPath, pngPath)
% PLOT_IV Generate publication-quality I-V characteristic plot.

% Figure: single-column width (~8.6 cm for PRL/PRB), 4:3 aspect
fig = figure('Units', 'centimeters', 'Position', [2 2 12 9], ...
    'PaperUnits', 'centimeters', 'PaperSize', [12 9], ...
    'PaperPosition', [0 0 12 9], 'Color', 'w', 'Visible', 'off');
ax = axes(fig);
hold(ax, 'on');

% Colors: professional palette
c_tol = [0.1 0.1 0.1];      % near-black for total
c_seq = [0.0 0.35 0.75];    % blue for sequential
c_cot = [0.85 0.15 0.15];   % red for cotunneling

% Plot: markers for total, solid lines for components
plot(ax, Vsd, I_tol, 'o', ...
    'Color', c_tol, 'MarkerSize', 2.5, 'MarkerFaceColor', c_tol, ...
    'DisplayName', '$I_{\mathrm{total}}$');
plot(ax, Vsd, I_seq, '-', ...
    'Color', c_seq, 'LineWidth', 1.4, ...
    'DisplayName', '$I_{\mathrm{seq}}$');
plot(ax, Vsd, I_cot, '--', ...
    'Color', c_cot, 'LineWidth', 1.4, ...
    'DisplayName', '$I_{\mathrm{cot}}$');

% Labels (LaTeX)
xlabel(ax, '$V_{\mathrm{sd}}$ (V)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel(ax, '$I$ (A)', 'Interpreter', 'latex', 'FontSize', 12);

% Title with key parameters
titleStr = sprintf('$N{=}%d,\\;\\lambda{=}%.1f,\\;T{=}%.1f$ K,\\;$\\hbar\\omega{=}%.0f$ meV', ...
    N, lambda, T, vmode*1e3);
title(ax, titleStr, 'Interpreter', 'latex', 'FontSize', 11);

% Legend
lgd = legend(ax, 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'northwest');
lgd.BoxFace.ColorType = 'truecoloralpha';
lgd.BoxFace.ColorData = uint8([255 255 255 220]');

% Axes styling
ax.FontSize = 10;
ax.TickLabelInterpreter = 'latex';
ax.LineWidth = 0.6;
box(ax, 'on');
grid(ax, 'on');
ax.GridLineStyle = '-';
ax.GridAlpha = 0.12;
ax.MinorGridLineStyle = 'none';

% Annotation: wall time + parameters
annoStr = sprintf('$\\Gamma_{L,R}/\\hbar\\omega = %.2f$, $\\eta = %.1f$, $V_g = %.1f$ V\\newline MATLAB, $t = %.1f$ s', ...
    alphaL, eta, Vg, wall_time);
annotation(fig, 'textbox', [0.48 0.12 0.48 0.16], ...
    'String', annoStr, 'Interpreter', 'latex', ...
    'FontSize', 8, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
    'FitBoxToText', 'on', 'BackgroundColor', [1 1 1 0.7]);

hold(ax, 'off');

% Save
try
    exportgraphics(fig, pdfPath, 'ContentType', 'vector', 'BackgroundColor', 'white');
    exportgraphics(fig, pngPath, 'Resolution', 300, 'BackgroundColor', 'white');
catch
    % Fallback for MATLAB < R2020a
    print(fig, pdfPath, '-dpdf', '-r300');
    print(fig, pngPath, '-dpng', '-r300');
end

close(fig);
fprintf('Plot saved.\n');

end
