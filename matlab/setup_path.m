function setup_path()
% SETUP_PATH Add all required directories to the MATLAB path
%   Run this once per session before using the Franck-Condon simulation.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, 'src'));
addpath(fullfile(thisDir, 'utils'));
addpath(fullfile(thisDir, 'test'));

fprintf('Franck-Condon simulation paths added.\n');
end
