%SETUP_PATHS  Add the shared helpers in utils/ to the MATLAB path.
%
%  Run once per MATLAB session from the repository root:
%       >> setup_paths
%
%  EEGLAB must also be on the path (run `eeglab` once, or addpath it).

repo_root = fileparts(mfilename('fullpath'));
addpath(fullfile(repo_root, 'utils'));

if exist('pop_loadset', 'file') ~= 2
    warning('setup_paths:noEEGLAB', ...
        'EEGLAB not found on the path. Install it and run `eeglab` before the exercises.');
end

fprintf('Signal Processing Exercises: utils/ added to path.\n');
clear repo_root
