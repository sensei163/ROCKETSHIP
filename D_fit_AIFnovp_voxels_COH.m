
%% 3. Fit to model No VP
function a = D_fit_AIFnovp_voxels()
%% a) Load files from C_fit_AIFnoCP.m
%

r2filter = 0; % Filter out all fits with r2 < r2filter
CPU      = 0; % Use only if you are doing multicore


% Load the data files

rootname  = 'new';
[gogo,PathName1,FilterIndex] = uigetfile(['/data/studies/' '/*AIF_no_vpFIT_ROI.mat'],'Choose R1 file');

% Load the data files
directory = PathName1
rootname  = strrep(gogo, '.nii', '');
load(fullfile(PathName1, gogo));


xdata{1}.numvoxels = numvoxels;

% %% b) voxel by voxel fitting
% 
% neuroecon = 1;
% tic
% if(neuroecon)
%     
%     %Schedule object, neuroecon
%     
%     sched = findResource('scheduler', 'configuration', 'NeuroEcon.local')
%     set(sched, 'SubmitArguments', '-l walltime=2:00:00 -m abe -M thomasn@caltech.edu')
%     
%     
%     warning off
%     
%     p = pwd;
%     n = '/home/thomasn/scripts/niftitools';
% 
%     
%     
%     
%     %j = createMatlabPoolJob(sched, 'configuration', 'NeuroEcon.local','FileDependencies', {fullfile(p, 'FXLStep1AIFhelper.m')})
%     j = createMatlabPoolJob(sched, 'configuration', 'local','PathDependencies', {p})
%     
%     set(j, 'MaximumNumberOfWorkers', CPU)
%     set(j, 'MinimumNumberOfWorkers', 1)
%     
%     createTask(j, @FXLfit_novpC, 1,{xdata, xdata{1}.numvoxels});
%     
%     
%     
%     submit(j);
%     waitForState(j);
%     
%     
%     results = getAllOutputArguments(j);
%  
%     destroy(j);
%     %getDebugLog(sched, j)
%     %pause
%     x = cell2mat(results);
% else
    
    x = FXLfit_novpC(xdata, xdata{1}.numvoxels);
%end


toc
%save(fullfile(fileparts(directory1), 'temptemp.mat'))
%% c.1 Check the fit
if(CHECK)
    check = CHECKnovp(x, xdata);
end


save(fullfile(fileparts(directory1), [rootname 'AIF_no_vpFIT_voxels.mat']), 'x', 'tumind','dynamname', 'directory1' , 'xdata')

%% d) Check if physiologically possible, if not, remove

% ve > 1
checkind = find(x(:,2) > 1);
x(checkind,:) = [];
tumind(checkind) = [];



% ve < 0
checkind = find(x(:,2) < 0);
x(checkind,:) = [];
tumind(checkind) = [];


% Ktrans > 3
checkind = find(x(:,1) > 3);
x(checkind,:) = [];
tumind(checkind) = [];

% Ktrans < 0
checkind = find(x(:,1) < 0);
x(checkind,:) = [];
tumind(checkind) = [];

%% d.2) Check R2 fit

checkind = find(x(:,4) < r2filter);
x(checkind,:)    = [];
tumind(checkind) = [];

Ktrans = x(:,1);
ve     = x(:,2);
r2     = x(:,3);


%% f) Make maps

KtransROI = zeros(size(currentimg));
veROI     = zeros(size(currentimg));
R2        = zeros(size(currentimg));

size(tumind)
size(Ktrans)

KtransROI(tumind) = Ktrans;
veROI(tumind)     = ve;
R2(tumind)        = r2;



%% g) Save file
res = [0 0.25 0.25 2];


%% h) Save image files

[discard actual] = fileparts(strrep(dynamname.fileprefix, '\', '/'));

save_nii(make_nii(KtransROI, res(2:4), [1 1 1]), fullfile(fileparts(directory1), [actual '_Ktrans_novp.nii']));
save_nii(make_nii(veROI, res(2:4), [1 1 1]), fullfile(fileparts(directory1),[actual '_ve_novp.nii']));
save_nii(make_nii(R2, res(2:4), [1 1 1]), fullfile(fileparts(directory1),[actual '_R2_novp.nii']));


a = 1;



