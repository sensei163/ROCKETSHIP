
%% 3. Fit to model No VP
function a = D_fit_FXR_voxels_COH(directory1, CPU, CHECK, r2filter)
%% a) Load files from B_fit_AIFnoCP.m
%
%r2filter = 0; % Filter out all fits with r2 < r2filter
%CPU      = 0; % Use only if you are doing multicore


% Load the data files

rootname  = 'parker';
% [gogo,PathName1,FilterIndex] = uigetfile(['/data/studies/' '/*AIF_FXR_ROI.mat'],'Choose R1 file');
% 
% % Load the data files
% directory = PathName1
% rootname  = strrep(gogo, '.nii', '');
% load(fullfile(PathName1, gogo));
load(directory1)


xdata{1}.numvoxels = numvoxels;

%% b) voxel by voxel fitting

neuroecon = 1;
tic
if(neuroecon)
    
    %Schedule object, neuroecon
    
    sched = findResource('scheduler', 'configuration', 'NeuroEcon.local')
    set(sched, 'SubmitArguments', '-l walltime=5:00:00 -m abe -M thomasn@caltech.edu')
    
    
    warning off
    
    p = pwd
    n = '/home/thomasn/scripts/niftitools';
    
    
    
    j = createMatlabPoolJob(sched, 'configuration', 'NeuroEcon.local','PathDependencies', {p});
     set(j, 'MinimumNumberOfWorkers', 1)
    set(j, 'MaximumNumberOfWorkers', CPU)
   
    
    createTask(j, @FXRfit, 1,{xdata, xdata{1}.numvoxels})

    submit(j);
    waitForState(j);

    results = getAllOutputArguments(j)
    j
    destroy(j);

    %pause
    x = cell2mat(results);
else
    
    x = FXRfit(xdata, xdata{1}.numvoxels);
end

toc
%save(fullfile(fileparts(directory1), 'temptemp.mat'))
%% c.1 Check the fit
if(CHECK)
    check = CHECKFXR(x, xdata);
end


save(fullfile(fileparts(directory1), [rootname 'AIF_FXR_voxels.mat']), 'x', 'tumind','dynamname', 'directory1' , 'xdata');

%% d) Check if physiologically possible, if not, remove

%% Error

checkind = find(x(:,1) < 0);
x(checkind,:) = [];
tumind(checkind) = [];

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
tau    = x(:,3);
r2     = x(:,4);


%% f) Make maps

KtransROI = zeros(size(currentimg));
veROI     = zeros(size(currentimg));
tauROI    = zeros(size(currentimg));
R2        = zeros(size(currentimg));

size(tumind)
size(Ktrans)

KtransROI(tumind) = Ktrans;
veROI(tumind)     = ve;
tauROI(tumind)    = tau;
R2(tumind)        = r2;



%% g) Save file
res = [0 0.25 0.25 2];


%% h) Save image files

[discard actual] = fileparts(strrep(dynamname.fileprefix, '\', '/'));

save_nii(make_nii(KtransROI, res(2:4), [1 1 1]), fullfile(fileparts(directory1), [actual '_Ktrans_FXR.nii']));
save_nii(make_nii(veROI, res(2:4), [1 1 1]), fullfile(fileparts(directory1),[actual '_ve_FXR.nii']));
save_nii(make_nii(tauROI, res(2:4), [1 1 1]), fullfile(fileparts(directory1),[actual '_tau_FXR.nii']));
save_nii(make_nii(R2, res(2:4), [1 1 1]), fullfile(fileparts(directory1),[actual '_R2_FXR.nii']));


a = 1;



