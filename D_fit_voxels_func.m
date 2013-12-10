% 3. Fit to model VP on a voxel by voxel basis

function results = D_fit_voxels_func(results_b_path,dce_model,time_smoothing,time_smoothing_window,xy_smooth_size,number_cpus,neuroecon)

% Toggle options
%************************
r2filter = 0;		% Filter out all fits with r2 < r2filter

% number_cpus = 4;	% Use only if you are doing multicore

close_pool = 0;		% Close matlabpool when done with processing

% base_directory = 'C:\Users\sbarnes\Documents\data\6 DCE Stroke\sb01_06nov13.mH1';

% neuroecon =0;

% time_smoothing = 2;	%0=none, 1=moving, 2=rlowess

% Select fitting model
% 'aif_vp' = tofts with vascular compartment
% 'aif' = tofts without vascular compartment
% 'fxr' = not implemented
% 'sauc' = not implemented 
% 'ss' = not implemented 
% 'fractal' = not implemented
% 'auc' = not implemented
% 'auc_rr' = not implemented
% dce_model = 'aif';

% End options
%************************

% a) Load the data files
load(results_b_path);

[PathName,base,ext] = fileparts(results_b_path);
rootname  = base;
rootname  = strrep(rootname, 'fitted_R1info', '');

xdata{1}.numvoxels = numvoxels;
disp(['Fitting data using the ' dce_model ' model']);
disp(['Staring DCE processing on ' num2str(numvoxels) ' voxels...']);

% Open pool if not open or improperly sized
if matlabpool('size')~= number_cpus
	if matlabpool('size')>0
		matlabpool close;
	end
	matlabpool('local', number_cpus);
end

% Save original data
xdata{1}.Ct_original = xdata{1}.Ct;

% a1) smoothing in image domain
if xy_smooth_size~=0
	% original_timepoint = NaN(size(currentimg));
	original_timepoint = zeros(size(currentimg));
	h = fspecial('gaussian', [xy_smooth_size xy_smooth_size],0.5);
	for i=1:size(xdata{1}.Ct,1)
		original_timepoint(tumind) = xdata{1}.Ct(i,:);
	% 	imshow(smooth_timepoint.*20000)
		smooth_timepoint = filter2(h, original_timepoint);
		xdata{1}.Ct(i,:) = smooth_timepoint(tumind)';
	end
end

% a2) smoothing in time domain
if ~strcmp(time_smoothing,'none')
	disp('Smoothing time domain');
	
	Ct_all = xdata{1}.Ct;
	p = ProgressBar(size(Ct_all,2));
	parfor i=1:size(Ct_all,2)
		Ct_smooth = Ct_all(:,i);
		if strcmp(time_smoothing,'moving')
			Ct_smooth = smooth(Ct_smooth,time_smoothing_window,'moving');
		elseif strcmp(time_smoothing,'rlowess')
			Ct_smooth = smooth(Ct_smooth,time_smoothing_window/size(Ct_smooth,1),'rlowess');
		else
			% no smoothing
		end

		Ct_all(:,i) = Ct_smooth;
		p.progress;
	end
	xdata{1}.Ct = Ct_all;
	p.stop;
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
        t = createTask(jj, @FXLfit_generic, 1,{xdata, numvoxels, dce_model,time_smoothing,time_smoothing_window});%numel(STARTEND(k,1):STARTEND(k,2))});
        set(t, 'CaptureCommandWindowOutput', true);
   
        submit(jj)
        waitForState(jj,'finished')
        jj
        results = getAllOutputArguments(jj)
        destroy(jj)
    
        clear jj
        fitting_results = cell2mat(results);
       % x(STARTEND(k,1):STARTEND(k,2),:) = cell2mat(results);  
else
	disp('Fitting Voxels');
    fitting_results = FXLfit_generic(xdata, numvoxels, dce_model);
end
processing_time = toc;
disp(['processing completed in ' datestr(processing_time/86400, 'HH:MM:SS') ' (hr:min:sec)']);
if close_pool
	matlabpool close;
end

% c) Save file
%************************
save(fullfile(PathName, [rootname dce_model '_fit_voxels.mat']), 'fitting_results', 'tumind','dynamname', 'PathName' , 'xdata','time_smoothing','time_smoothing_window', 'dce_model')
results = fullfile(PathName, [rootname dce_model '_fit_voxels.mat']);

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

if strcmp(dce_model, 'aif')
	KtransROI = zeros(size(currentimg));
	veROI     = zeros(size(currentimg));
	residual  = zeros(size(currentimg));
	ci_95_low_ktrans	= zeros(size(currentimg));
	ci_95_high_ktrans	= zeros(size(currentimg));
	ci_95_low_ve		= zeros(size(currentimg));
	ci_95_high_ve		= zeros(size(currentimg));	

	KtransROI(tumind) = fitting_results(:,1);
	veROI(tumind)     = fitting_results(:,2);
	residual(tumind)  = fitting_results(:,4);
	ci_95_low_ktrans(tumind)	= fitting_results(:,5);
	ci_95_high_ktrans(tumind)	= fitting_results(:,6);
	ci_95_low_ve(tumind)		= fitting_results(:,7);
	ci_95_high_ve(tumind)		= fitting_results(:,8);
	

	save_nii(make_nii(KtransROI, res(2:4), [1 1 1]), fullfile(PathName, [rootname dce_model '_Ktrans.nii']) );
	save_nii(make_nii(veROI, res(2:4), [1 1 1]), fullfile(PathName, [rootname dce_model '_ve.nii']));
	save_nii(make_nii(residual, res(2:4), [1 1 1]), fullfile(PathName, [rootname dce_model '_residual.nii']));
	save_nii(make_nii(ci_95_low_ktrans, res(2:4), [1 1 1]), fullfile(PathName, [rootname dce_model '_ktrans_ci_low.nii']));
	save_nii(make_nii(ci_95_high_ktrans, res(2:4), [1 1 1]), fullfile(PathName, [rootname dce_model '_ktrans_ci_high.nii']));
	save_nii(make_nii(ci_95_low_ve, res(2:4), [1 1 1]), fullfile(PathName, [rootname dce_model '_ve_ci_low.nii']));
	save_nii(make_nii(ci_95_high_ve, res(2:4), [1 1 1]), fullfile(PathName, [rootname dce_model '_ve_ci_high.nii']));

elseif strcmp(dce_model, 'aif_vp')
	KtransROI = zeros(size(currentimg));
	veROI     = zeros(size(currentimg));
	vpROI     = zeros(size(currentimg));
	residual  = zeros(size(currentimg));
	ci_95_low_ktrans	= zeros(size(currentimg));
	ci_95_high_ktrans	= zeros(size(currentimg));
	ci_95_low_ve		= zeros(size(currentimg));
	ci_95_high_ve		= zeros(size(currentimg));	
	ci_95_low_vp		= zeros(size(currentimg));
	ci_95_high_vp		= zeros(size(currentimg));	
	
	KtransROI(tumind) = fitting_results(:,1);
	veROI(tumind)     = fitting_results(:,2);
	vpROI(tumind)     = fitting_results(:,3);
	residual(tumind)  = fitting_results(:,4);
	ci_95_low_ktrans(tumind)	= fitting_results(:,5);
	ci_95_high_ktrans(tumind)	= fitting_results(:,6);
	ci_95_low_ve(tumind)		= fitting_results(:,7);
	ci_95_high_ve(tumind)		= fitting_results(:,8);
	ci_95_low_vp(tumind)		= fitting_results(:,9);
	ci_95_high_vp(tumind)		= fitting_results(:,10);
	
	save_nii(make_nii(KtransROI, res(2:4), [1 1 1]), fullfile(PathName, [rootname dce_model '_Ktrans.nii']));
	save_nii(make_nii(veROI, res(2:4), [1 1 1]), fullfile(PathName, [rootname dce_model '_ve.nii']));
	save_nii(make_nii(vpROI, res(2:4), [1 1 1]), fullfile(PathName, [rootname dce_model '_vp.nii']));
	save_nii(make_nii(residual, res(2:4), [1 1 1]), fullfile(PathName, [rootname dce_model '_residual.nii']));
	save_nii(make_nii(ci_95_low_ktrans, res(2:4), [1 1 1]), fullfile(PathName, [rootname dce_model '_ktrans_ci_low.nii']));
	save_nii(make_nii(ci_95_high_ktrans, res(2:4), [1 1 1]), fullfile(PathName, [rootname dce_model '_ktrans_ci_high.nii']));
	save_nii(make_nii(ci_95_low_ve, res(2:4), [1 1 1]), fullfile(PathName, [rootname dce_model '_ve_ci_low.nii']));
	save_nii(make_nii(ci_95_high_ve, res(2:4), [1 1 1]), fullfile(PathName, [rootname dce_model '_ve_ci_high.nii']));
	save_nii(make_nii(ci_95_low_vp, res(2:4), [1 1 1]), fullfile(PathName, [rootname dce_model '_vp_ci_low.nii']));
	save_nii(make_nii(ci_95_high_vp, res(2:4), [1 1 1]), fullfile(PathName, [rootname dce_model '_vp_ci_high.nii']));
end

disp('Finished D');




