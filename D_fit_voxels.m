% 3. Fit to model VP on a voxel by voxel basis

% Load files from C_fit_AIFnoCP.m

% Toggle options
%************************
r2filter = 0;		% Filter out all fits with r2 < r2filter

number_cpus = 4;	% Use only if you are doing multicore

close_pool = 0;		% Close matlabpool when done with processing

base_directory = 'C:\Users\sbarnes\Documents\data\6 DCE Stroke\sb01_06nov13.mH1';

neuroecon =0;

smooth_model = 2;	%0=none, 1=moving, 2=lowess

% Select fitting model
% 'aif_vp' = tofts with vascular compartment
% 'aif' = tofts without vascular compartment
% 'fxr' = not implemented
% 'sauc' = not implemented 
% 'ss' = not implemented 
% 'fractal' = not implemented
% 'auc' = not implemented
% 'auc_rr' = not implemented
model = 'aif';

% End options
%************************

% a) Load the data files
[gogo,PathName,FilterIndex] = uigetfile([base_directory '/*AIF_with_vpFIT_ROI.mat'],'Choose R1 file');
if gogo==0
	disp('User selected cancel');
	return
end
directory = PathName;
rootname  = strrep(gogo, '.nii', '');
load(fullfile(PathName, gogo));

xdata{1}.numvoxels = numvoxels;
disp(['Fitting data using the ' model ' model']);
disp(['Staring DCE processing on ' num2str(numvoxels) ' voxels...']);


% a1) smoothing in image domain
% original_timepoint = NaN(size(currentimg));
original_timepoint = zeros(size(currentimg));
h = fspecial('average', [2 2]);
for i=1:size(xdata{1}.Ct,1)
	original_timepoint(tumind) = xdata{1}.Ct(i,:);
% 	imshow(smooth_timepoint.*20000)
	smooth_timepoint = filter2(h, original_timepoint);
	xdata{1}.Ct(i,:) = smooth_timepoint(tumind)';
end

% b) voxel by voxel fitting
%************************
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
%         STARTEND(k,:)
%         %We only feed the workers only the voxels that they can handle
%         
%         xdata{1}.Ct = wholeCt(:,STARTEND(k,1):STARTEND(k,2));
        
        %Schedule object, neuroecon
        t = createTask(jj, @FXLfit_generic, 1,{xdata, numvoxels, model,smooth_model});%numel(STARTEND(k,1):STARTEND(k,2))});
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
	% Open pool if not open or improperly sized
	if matlabpool('size')~= number_cpus
		matlabpool close;
		matlabpool('local', number_cpus);
	end
    x = FXLfit_generic(xdata, numvoxels, model,smooth_model);
	if close_pool
		matlabpool close;
	end
end
processing_time = toc;
disp(['processing completed in ' datestr(processing_time/86400, 'HH:MM:SS') ' (hr:min:sec)']);


% c) Save file
%************************
save(fullfile(fileparts(directory), [rootname model '_FIT_voxels.mat']), 'x', 'tumind','dynamname', 'directory' , 'xdata')

% d) Check if physiologically possible, if not, remove
%************************
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


% f) Make maps and Save image files
%************************
[discard, actual] = fileparts(strrep(dynamname.fileprefix, '\', '/'));
res = [0 0.25 0.25 2];

if strcmp(model, 'aif')
	KtransROI = zeros(size(currentimg));
	veROI     = zeros(size(currentimg));
	residual  = zeros(size(currentimg));

	KtransROI(tumind) = x(:,1);
	veROI(tumind)     = x(:,2);
	residual(tumind)  = x(:,4);

	save_nii(make_nii(KtransROI, res(2:4), [1 1 1]), fullfile(fileparts(directory), [actual '_' model '_Ktrans.nii']));
	save_nii(make_nii(veROI, res(2:4), [1 1 1]), fullfile(fileparts(directory),[actual '_' model '_ve.nii']));
	save_nii(make_nii(residual, res(2:4), [1 1 1]), fullfile(fileparts(directory),[actual '_' model '_residual.nii']));
elseif strcmp(model, 'aif_vp')
	KtransROI = zeros(size(currentimg));
	veROI     = zeros(size(currentimg));
	vpROI     = zeros(size(currentimg));
	residual  = zeros(size(currentimg));

	KtransROI(tumind) = x(:,1);
	veROI(tumind)     = x(:,2);
	vpROI(tumind)     = x(:,3);
	residual(tumind)  = x(:,4);

	save_nii(make_nii(KtransROI, res(2:4), [1 1 1]), fullfile(fileparts(directory), [actual '_' model '_Ktrans.nii']));
	save_nii(make_nii(veROI, res(2:4), [1 1 1]), fullfile(fileparts(directory),[actual '_' model '_ve.nii']));
	save_nii(make_nii(vpROI, res(2:4), [1 1 1]), fullfile(fileparts(directory),[actual '_' model '_vp.nii']));
	save_nii(make_nii(residual, res(2:4), [1 1 1]), fullfile(fileparts(directory),[actual '_' model '_residual.nii']));
end






