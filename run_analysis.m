% Find and add subpaths 
 mfilepath=fileparts(which('run_dce'));
 addpath(fullfile(mfilepath,'dce'));
 addpath(fullfile(mfilepath,'external_programs'));
 addpath(fullfile(mfilepath,'external_programs/niftitools'));
 addpath(fullfile(mfilepath,'parametric_scripts'));
% Run program
fitting_analysis