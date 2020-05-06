% Find and add subpaths 
 mfilepath=fileparts(which('run_dsc'));
 addpath(fullfile(mfilepath,'dsc'));
 addpath(fullfile(mfilepath,'external_programs'));
 addpath(fullfile(mfilepath,'external_programs/niftitools'));
 addpath(fullfile(mfilepath,'parametric_scripts'));
% Run program
 DSC_gui
