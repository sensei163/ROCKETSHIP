% 3. Fit to model VP on a voxel by voxel basis

% function a = D_fit_AIFwithvp_voxels_COH(directory1, CPU, CHECK, r2filter)
% 1. Load files from C_fit_AIFnoCP.m

r2filter = 0; % Filter out all fits with r2 < r2filter
number_cpus = 4; % Use only if you are doing multicore
% CPU

% Load the data files

rootname  = 'parker';
[gogo,PathName1,FilterIndex] = uigetfile(['C:\Users\sbarnes\Documents\data\6 DCE Stroke\sb01_06nov13.mH1' '/*AIF_with_vpFIT_ROI.mat'],'Choose R1 file');

% Load the data files
directory = PathName1
rootname  = strrep(gogo, '.nii', '');
load(fullfile(PathName1, gogo));
% load(directory1)





xdata{1}.numvoxels = numvoxels;

% totale = ceil(numvoxels/CPU);
% %Split into rows by CPU.
% for i = 1:totale
%     STARTEND(i,1) = (i-1)*CPU+1;
%     STARTEND(i,2) = min(numvoxels, i*CPU);
% end

%wholeCt = xdata{1}.Ct;

%%% b) voxel by voxel fitting

%Example scripts for scheduler.

neuroecon =0;
tic
if(neuroecon)
    warning off
    
    p = pwd
    n = '/home/thomasn/scripts/niftitools';
   % for k = 1:totale 
        sched = findResource('scheduler', 'configuration', 'NeuroEcon.local');
        set(sched, 'SubmitArguments', '-l walltime=5:00:00 -m abe -M thomasn@caltech.edu')
        
        jj = createMatlabPoolJob(sched, 'PathDependencies', {p});
        
        set(jj, 'MaximumNumberOfWorkers', number_cpus)
        set(jj, 'MinimumNumberOfWorkers', number_cpus)      
%         
%         STARTEND(k,:)
%         %We only feed the workers only the voxels that they can handle
%         
%         xdata{1}.Ct = wholeCt(:,STARTEND(k,1):STARTEND(k,2));
        
        %Schedule object, neuroecon
        
        t = createTask(jj, @FXLfit_withvpC, 1,{xdata, numvoxels});%numel(STARTEND(k,1):STARTEND(k,2))});
        set(t, 'CaptureCommandWindowOutput', true);
   
        submit(jj)
        
 
        waitForState(jj,'finished')
        jj
        results = getAllOutputArguments(jj)
        destroy(jj)
    
        clear jj
        x = cell2mat(results);
       % x(STARTEND(k,1):STARTEND(k,2),:) = cell2mat(results);  
else
	matlabpool('local', number_cpus);
    x = FXLfit_withvpC(xdata, numvoxels);
    size(x)
	matlabpool close;
end
toc

%% g) Save file

save(fullfile(fileparts(directory), [rootname 'AIF_with_vpFIT_voxels.mat']), 'x', 'tumind','dynamname', 'directory' , 'xdata')

%% d) Check if physiologically possible, if not, remove

%% Error

% checkind = find(x(:,1) < 0);
% x(checkind,:) = [];
% tumind(checkind) = [];

% ve > 1
% checkind = find(x(:,2) > 1);
% x(checkind,:) = [];
% tumind(checkind) = [];
% 
% % ve < 0
% checkind = find(x(:,2) < 0);
% x(checkind,:) = [];
% tumind(checkind) = [];
% 
% % vp > 1
% checkind = find(x(:,3) > 1);
% x(checkind,:) = [];
% tumind(checkind) = [];
% 
% % vp < 0
% checkind = find(x(:,3) < 0);
% x(checkind,:) = [];
% tumind(checkind) = [];
% 
% % Ktrans > 3
% checkind = find(x(:,1) > 3);
% x(checkind,:) = [];
% tumind(checkind) = [];
% 
% % Ktrans < 0
% checkind = find(x(:,1) < 0);
% x(checkind,:) = [];
% tumind(checkind) = [];
% 
% %% d.2) Check R2 fit
% 
% checkind = find(x(:,4) < r2filter);
% x(checkind,:)    = [];
% tumind(checkind) = [];

Ktrans = x(:,1);
ve     = x(:,2);
vp     = x(:,3);
r2     = x(:,4);


%% f) Make maps

KtransROI = zeros(size(currentimg));
veROI     = zeros(size(currentimg));
vpROI     = zeros(size(currentimg));
R2        = zeros(size(currentimg));


KtransROI(tumind) = Ktrans;
veROI(tumind)     = ve;
vpROI(tumind)     = vp;
R2(tumind)        = r2;

res = [0 0.25 0.25 2];

%% h) Save image files

[discard actual] = fileparts(strrep(dynamname.fileprefix, '\', '/'));

save_nii(make_nii(KtransROI, res(2:4), [1 1 1]), fullfile(fileparts(directory), [actual '_Ktrans_withvp.nii']));
save_nii(make_nii(veROI, res(2:4), [1 1 1]), fullfile(fileparts(directory),[actual '_ve_withvp.nii']));
save_nii(make_nii(vpROI, res(2:4), [1 1 1]), fullfile(fileparts(directory),[actual '_vp_withvp.nii']));
save_nii(make_nii(R2, res(2:4), [1 1 1]), fullfile(fileparts(directory),[actual '_R2_novp.nii']));

a = 1;





