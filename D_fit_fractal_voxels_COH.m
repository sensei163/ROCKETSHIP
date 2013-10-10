%% 3. Fit to model VP on a voxel by voxel basis

function xdata = D_fit_fractal_voxels_COH(directory1, CPU, CHECK, r2filter, timestamp)
%% 1. Load files from C_fit_AIFnoCP.m

% r2filter = 0; % Filter out all fits with r2 < r2filter
% CPU      = 0; % Use only if you are doing multicore
CPU

% Load the data files

rootname  = 'AUC';
% [gogo,PathName1,FilterIndex] = uigetfile(['/data/studies/' '/*AIF_with_vpFIT_ROI.mat'],'Choose R1 file');
%
% % Load the data files
% directory = PathName1
% rootname  = strrep(gogo, '.nii', '');
% load(fullfile(PathName1, gogo));
load(directory1)



%% f) Make maps

AUCROI    = zeros(size(currentimg));


AUCROI(tumind) = NAUC;

xdata{1}.AUC = AUC;
xdata{1}.NAUC= NAUC;
xdata{1}.AUCCp=AUCCp;


res = [0 0.25 0.25 2];

%% h) Save image files

[discard actual] = fileparts(strrep(dynamname.fileprefix, '\', '/'));

save_nii(make_nii(AUCROI, res(2:4), [1 1 1]), fullfile(fileparts(directory1), [actual '_NAUC_min' num2str(timestamp) '.nii']));

a = 1;





