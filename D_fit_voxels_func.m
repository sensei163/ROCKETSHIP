function [results, batch] = D_fit_voxels_func(results_b_path,dce_model,time_smoothing,time_smoothing_window,xy_smooth_size,number_cpus,roi_list,fit_voxels,neuroecon, batch, outputft)

% D_fit_voxels_func - Fit DCE curve to various models on a voxel by voxel
% or ROI basis
%
% Inputs:
%  results_b_path     - *.mat Results from part B
%  dce_model          - Select fitting model
%                       'aif_vp' = tofts with vascular compartment
%                       'aif' = tofts without vascular compartment
%                       'fxr' = "shutter speed"
%                       'sauc' = not implemented
%                       'ss' = not implemented
%                       'fractal' = not implemented
%                       'auc' = not implemented
%                       'auc_rr' = not implemented
%  time_smoothing     - type of time smoothing
%                       'none' = no smoothing
%                       'moving' = moving average
%                       'rlowess' = robust local regression
%  time_smoothing_window - size of time smoothing window (time points)
%  xy_smooth_size     -	sigma of the Gaussian low pass smooth function
%  number_cpus        - number of cpu cores for parallel processing
%  roi_list           - paths to ROIs that specify homogenous regions to
%                       calculate a single DCE fit for. Values >0
%                       considered in the ROI, values <=0 considered
%                       outside the ROI
%  fit_voxels         - perform DCE fit on individual voxels
%  neuroecon          - perform processing on neuroecon server
%  batch              - prep for batch procesing only
%  outputft           - output filetype
%
% The script loads the data arrays generated from B_AIF_fitting_func().
% Then it will fit a DCE curve according to various models
%
% Requires:
% D_fit_voxels_func.m
% DataHash.m
% FXLfit_generic.m
% ProgressBar.m
% FXLStep1AIFhelper.m
% FXLStep1AIFhelper_vp.m
% fxr_helper.m
% parse_preference_file.m
% niftitools
%
% Samuel Barnes
% Caltech
% December 2013
% Updated Thomas Ng Jan 2014


% Toggle options
%************************
r2filter = 0;		% Filter out all fits with r2 < r2filter
close_pool = 0;		% Close matlabpool when done with processing
% End options
%************************

% Setup DCE model string

dce_model_string = '';

if dce_model.aif
    dce_model_string = [dce_model_string 'Tofts'];
end
if dce_model.aif_vp
    dce_model_string = [dce_model_string ', Tofts w/ Vp'];
end
if dce_model.fxr
    dce_model_string = [dce_model_string ', FXR'];
end

if dce_model.fractal
    dce_model_string = [dce_model_string ', Fractal metric'];
end
if dce_model.auc
    dce_model_string = [dce_model_string ', Area under curve'];
end

% a) Load the data files
load(results_b_path);
results_a_path = Bdata.results_a_path;
load(results_a_path);

rootname    = Adata.rootname;
xdata       = Bdata.xdata;
R1tTOI      = Adata.R1tTOI;
start_time  = Bdata.start_time;
end_time    = Bdata.end_time;
T1TUM       = Adata.T1TUM;
tumind      = Adata.tumind;
dynam_name  = Adata.dynam_name;
numvoxels   = Bdata.numvoxels;
currentimg  = Adata.currentCT;
res         = Adata.res;
relaxivity  = Adata.relaxivity;
hdr         = Adata.hdr;
sliceloc    = Adata.sliceloc;

if start_time == 0
    start_time = 1;
end
if end_time == 0
    end_time = size(R1tTOI,1);
end

% Load and prep raw signal files if required

if dce_model.auc
    Sss                         = Adata.Sss;
    Ssstum                      = Adata.Ssstum;
    Stlv                        = Adata.Stlv;
    Stlv                        = Stlv(start_time:end_time,:);% This was not done in RunB as per others.
    Sttum                       = Adata.Sttum;
    Sttum                       = Sttum(start_time:end_time,:);% This was not done in RunB as per others.
    xdata{1}.start_injection    = Bdata.start_injection;
    xdata{1}.end_injection      = Bdata.end_injection;
    xdata{1}.Stlv               = Stlv;
    xdata{1}.Ssstum             = Ssstum;
    xdata{1}.Sttum              = Sttum;
    xdata{1}.Sss                = Sss;
    
end



% Load the variables as needed.


% update output path to be same as location of input
[PathName,~,~] = fileparts(results_a_path);

% Log input results
log_path = fullfile(PathName, ['D_' rootname 'dce_fit_voxels.log']);
if exist(log_path, 'file')==2
    delete(log_path);
end
diary(log_path);
fprintf('************** User Input **************\n\n');
disp('User selected part B results: ');
disp(results_b_path);
Opt.Input = 'file';
b_md5 = DataHash(results_b_path, Opt);
fprintf('File MD5 hash: %s\n\n', b_md5)
disp('User selected dce model: ');
fprintf('%s\n\n',dce_model_string);
disp('User selected time smoothing model: ');
fprintf('%s\n\n',time_smoothing);
disp('User selected time smoothing window size: ');
disp(time_smoothing_window);
disp('User selected XY smooth size (sigma)');
disp(xy_smooth_size);
disp('User selected number of CPU cores');
disp(number_cpus);
disp('User selected ROI list');
[nrows,ncols]= size(roi_list);
for row=1:nrows
    disp(roi_list{row,:})
end
disp(' ');
disp('User selected fit individual voxels');
disp(fit_voxels);
disp('User selected use neuroecon');
disp(neuroecon);
fprintf('************** End User Input **************\n\n\n');

disp('Starting Part D - Fitting Voxels/ROIs')
disp(datestr(now))
disp(' ');
tic

% Start processing
xdata{1}.numvoxels = numvoxels;
%disp(['Fitting data using the ' 'dce' ' model']);

% Open pool if not open or improperly sized
if matlabpool('size')~= number_cpus
    % Do not launch pool with diary on, locks the log file
    diary off;
    if matlabpool('size')>0
        matlabpool close;
    end
    matlabpool('local', number_cpus);
    diary on;
end

% Substitute R1 data for concentration data in curve to fit
% FXR model fits to the R1 data directly, not concentrations
if dce_model.fxr
    xdata{1}.Ct = R1tTOI(start_time:end_time,:);
end

% Save original data
xdata{1}.Ct_original = xdata{1}.Ct;

% Perpare any ROIs
number_rois = 0;
if ~isempty(roi_list) && ~strcmp('No Files',cell2mat(roi_list(1)))
    %Sanitize list
    for m=size(roi_list,1):-1:1
        testfile=cell2mat(roi_list(m));
        if ~exist(testfile, 'file')
            % File does not exist.
            warning( 'File does not exist' );
            disp(testfile);
            return;
        end
        for n=(m-1):-1:1
            comparefile=roi_list(n);
            if strcmp(testfile,comparefile)
                disp( 'Removing duplicates' );
                disp(comparefile);
                roi_list(m)=[];
            end
        end
    end
    number_rois = size(roi_list,1);
    roi_name = [];
    %After sanitizing make sure we have some left
    if number_rois~=0
        [~, roi_name, roi_ext] = arrayfun(@(x) fileparts(x{:}), roi_list, 'UniformOutput', false);
        
        %Load ROI, find the selected voxels
        for r=number_rois:-1:1
            single_file=cell2mat(roi_list(r));
            
            single_roi = load_nii(single_file);
            single_roi = double(single_roi.img);
            roi_index{r}= find(single_roi > 0);
        end
        
        original_t1 = zeros(size(currentimg));
        roi_r1 = zeros(number_rois,1);
        original_timepoint = zeros(size(currentimg));
        roi_series = zeros(size(xdata{1}.Ct,1),number_rois);
        for t=1:size(xdata{1}.Ct,1)
            if dce_model.fxr
                original_timepoint(tumind) = R1tTOI(t,:);
            else
                original_timepoint(tumind) = xdata{1}.Ct(t,:);
                
                if dce_model.auc
                    original_timepoint_signal(tumind) = xdata{1}.Sttum(t,:);
                end
            end
            % 		original_timepoint(roi_index{r}) = 1e-4;
            % 		imshow(original_timepoint.*20000);
            %Average ROI voxels, insert into time series
            for r=number_rois:-1:1
                roi_series(t,r) = mean(original_timepoint(roi_index{r}));
                
                if dce_model.auc
                    roi_series_signal(t,r) = mean(original_timepoint_signal(roi_index{r}));
                end
            end
        end
        % For FXR
        original_t1(tumind) = T1TUM;
        for r=number_rois:-1:1
            roi_r1(r) = 1./mean(original_t1(roi_index{r}));
        end
        %make backup
        roi_series_original = roi_series;
    end
end

if ~fit_voxels && number_rois==0
    print('nothing to fit, select an ROI file or check "fit voxels"');
    return;
end

% for r=number_rois:-1:1
% 	figure;
% 	plot(xdata{1}.timer,roi_series(:,r));
% end
% return;

% a1) smoothing in image domain
if xy_smooth_size~=0 && fit_voxels
    % original_timepoint = NaN(size(currentimg));
    original_timepoint = zeros(size(currentimg));
    % make size 3*sigma rounded to nearest odd
    xy_smooth_size_odd = 2.*round((xy_smooth_size*3+1)/2)-1;
    if xy_smooth_size_odd<3
        xy_smooth_size_odd = 3;
    end
    h = fspecial('gaussian', [xy_smooth_size_odd xy_smooth_size_odd],xy_smooth_size);
    for i=1:size(xdata{1}.Ct,1)
        original_timepoint(tumind) = xdata{1}.Ct(i,:);
        % 	imshow(smooth_timepoint.*20000)
        smooth_timepoint = filter2(h, original_timepoint);
        xdata{1}.Ct(i,:) = smooth_timepoint(tumind)';
        
        if dce_model.auc
            original_timepoint_signal(tumind) = xdata{1}.Sttum(i,:);
            smooth_timepoint_signal = filter2(h, original_timepoint_signal);
            xdata{1}.Sttum(i,:) = smooth_timepoint_signal(tumind);
        end
        
        
    end
end

% a2) smoothing in time domain
if ~strcmp(time_smoothing,'none')
    disp('Smoothing time domain');
    
    if fit_voxels
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
        
        if dce_model.auc
            disp('Smoothing time domain for raw signals');
            Sttum_all = xdata{1}.Sttum;
            p = ProgressBar(size(Ct_all,2));
            parfor i=1:size(Ct_all,2)
                Ct_smooth = Sttum_all(:,i);
                if strcmp(time_smoothing,'moving')
                    Ct_smooth = smooth(Ct_smooth,time_smoothing_window,'moving');
                elseif strcmp(time_smoothing,'rlowess')
                    Ct_smooth = smooth(Ct_smooth,time_smoothing_window/size(Ct_smooth,1),'rlowess');
                else
                    % no smoothing
                end
                
                Sttum_all(:,i) = Ct_smooth;
                p.progress;
            end
            
            xdata{1}.Sttum = Sttum_all;
            p.stop;
            
        end
        
    end
    
    for r=1:number_rois
        roi_series(:,r);
        roi_smooth = roi_series(:,r);
        if strcmp(time_smoothing,'moving')
            roi_smooth = smooth(roi_smooth,time_smoothing_window,'moving');
        elseif strcmp(time_smoothing,'rlowess')
            roi_smooth = smooth(roi_smooth,time_smoothing_window/size(roi_smooth,1),'rlowess');
        else
            % no smoothing
        end
        
        roi_series(:,r) = roi_smooth;
        
        if dce_model.auc
            roi_smooth = roi_series_signal(:,r);
            if strcmp(time_smoothing,'moving')
                roi_smooth = smooth(roi_smooth,time_smoothing_window,'moving');
            elseif strcmp(time_smoothing,'rlowess')
                roi_smooth = smooth(roi_smooth,time_smoothing_window/size(roi_smooth,1),'rlowess');
            else
                % no smoothing
            end
            
            roi_series_signal(:,r) = roi_smooth;  
        end
    end
end

% a.b) Prep for batch.
    if dce_model.fxr

        Ddatabatch.T1TUM = T1TUM;
        Ddatabatch.relaxivity = relaxivity;
    end
      
    if number_rois
                Ddatabatch.roi_r1 = roi_r1;
        Ddatabatch.roi_series_original = roi_series_original;
        Ddatabatch.roi_series  = roi_series;
        
        if dce_model.auc
            Ddatabatch.roi_series_signal = roi_series_signal;
        end
    end
    
    Ddatabatch.number_cpus = number_cpus;
    Ddatabatch.xdata       = xdata;
    Ddatabatch.numvoxels   = numvoxels;
    Ddatabatch.dce_model   = dce_model;
    Ddatabatch.number_rois = number_rois;
    Ddatabatch.roi_name    = roi_name;
    Ddatabatch.roi_list    = roi_list;
    Ddatabatch.neuroecon   = neuroecon;
    Ddatabatch.close_pool  = close_pool;
    Ddatabatch.rootname    = rootname;
    Ddatabatch.hdr         = hdr;
    Ddatabatch.outputft    = outputft;
    Ddatabatch.sliceloc    = sliceloc;
    Ddatabatch.res         = res;
 
    Ddatabatch.fit_voxels  = fit_voxels;
    Ddatabatch.timind      = tumind;
    Ddatabatch.dynam_name  = dynam_name;
    Ddatabatch.PathName = PathName;
    Ddatabatch.time_smoothing = time_smoothing;
    Ddatabatch.smoothing_window =time_smoothing_window;
    Ddatabatch.currentimg = currentimg;
    Ddatabatch.results_a_path = results_a_path;
    Ddatabatch.results_b_path = results_b_path;
    
    if batch
        save(fullfile(PathName, ['D_' rootname 'prep' '_fit_voxels.mat']),  'Ddatabatch')
        results = fullfile(PathName, ['D_' rootname 'prep' '_fit_voxels.mat']);
        Opt.Input = 'file';
        mat_md5 = DataHash(results, Opt);
        disp(' ')
        disp('MAT results saved to: ')
        disp(results)
        disp(['File MD5 hash: ' mat_md5])
        disp(' ');
        disp('Prepped D for batch');
        disp(datestr(now))
        toc
        diary off;

        return;
    end
    
    % b) voxel by voxel fitting
    %************************
    % tic
    
    results = D_fit_voxels_batch_func(Ddatabatch);
    
    % Setup the list of models to use
    
%     cur_dce_model_list = {};
%     
% if dce_model.aif
%     cur_dce_model_list{end+1} = 'aif';
% end
% if dce_model.aif_vp
%     cur_dce_model_list{end+1} = 'aif_vp';
% end
% if dce_model.fxr
%     cur_dce_model_list{end+1} = 'fxr';
% end
% 
% if dce_model.fractal
%     cur_dce_model_list{end+1} = 'fractal';
% end
% if dce_model.auc
%     cur_dce_model_list{end+1} = 'auc';
% end
% 
%     % Run a for loop for every model.
%  for i = 1:numel(cur_dce_model_list)  
%      fit_voxels = 1; %% DEBUG
%      clear fitting_results;
%      cur_dce_model = cur_dce_model_list{i};
%      disp('  ');
%      disp(['Begin making maps for ' cur_dce_model '...']);
%      
%     if(neuroecon)
%         if strcmp(cur_dce_model, 'fxr')
%             xdata{1}.R1o = 1./T1TUM;
%             xdata{1}.R1i = 1./T1TUM;
%             xdata{1}.relaxivity = relaxivity;
%         end
%         
%         fitting_results = run_neuroecon_job(number_cpus, xdata, numvoxels, cur_dce_model);
%         % x(STARTEND(k,1):STARTEND(k,2),:) = cell2mat(results);
%     else
%         if number_rois~=0
%             disp(['Starting fitting for ' num2str(number_rois) ' ROIs...']);
%             
%             roi_data{1}.Cp = xdata{1}.Cp;
%             roi_data{1}.timer = xdata{1}.timer;
%             roi_data{1}.Ct = roi_series;
%             
%             if strcmp(cur_dce_model, 'fxr')
%                 roi_data{1}.R1o = roi_r1;
%                 roi_data{1}.R1i = roi_r1;
%                 roi_data{1}.relaxivity = relaxivity;
%             end
%             
%             if strcmp(cur_dce_model, 'auc')
%                 roi_data{1}.Sttum = roi_series_signal;
%             end
%             
%             roi_results = FXLfit_generic(roi_data, number_rois, cur_dce_model);
%             
%             disp('ROI fitting done')
%             disp(' ')    
%             
%         end
%         if fit_voxels
%             disp(['Starting fitting for ' num2str(numvoxels) ' voxels...']);
%             
%             if strcmp(cur_dce_model, 'fxr')
%                 xdata{1}.R1o = 1./T1TUM;
%                 xdata{1}.R1i = 1./T1TUM;
%                 xdata{1}.relaxivity = relaxivity;
%             end
%             
%             fitting_results = FXLfit_generic(xdata, numvoxels, cur_dce_model);
%             
%             disp('Voxel fitting done')
%         end
%     end
%     % processing_time = toc;
%     % disp(['processing completed in ' datestr(processing_time/86400, 'HH:MM:SS') ' (hr:min:sec)']);
%     if close_pool
%         matlabpool close;
%     end
%     
%     % c) Save file
%     %************************
%     fit_data.fit_voxels = fit_voxels;
%     fit_data.tumind = tumind;
%     fit_data.dynam_name = dynam_name;
%     fit_data.PathName = PathName;
%     fit_data.time_smoothing = time_smoothing;
%     fit_data.time_smoothing_window =time_smoothing_window;
%     fit_data.model_name =dce_model;
%     fit_data.number_rois = number_rois;
%     xdata{1}.dimensions = size(currentimg);
%     
%     if number_rois~=0
%         xdata{1}.roi_series = roi_series;
%         xdata{1}.roi_series_original = roi_series_original;
%         fit_data.roi_results = roi_results;
%         fit_data.roi_name = roi_name;
%         
%         if strcmp(cur_dce_model, 'fxr')
%             xdata{1}.roi_r1 = roi_r1;
%             xdata{1}.relaxivity = relaxivity;
%         end
%         
%         if strcmp(cur_dce_model, 'auc')
%             xdata{1}.roi_series_signal = roi_series_signal;
%         end
% 
%     end
%     if fit_voxels
%         fit_data.fitting_results  = fitting_results;
%     end
%     
%     Ddata.xdata = xdata;
%     Ddata.fit_data = fit_data;
%     Ddata.results_a_path = results_a_path;
%     Ddata.results_b_path = results_b_path;
%     
%     save(fullfile(PathName, ['D_' rootname cur_dce_model '_fit_voxels.mat']),  'Ddata')
%     results = fullfile(PathName, ['D_' rootname cur_dce_model '_fit_voxels.mat']);
%     Opt.Input = 'file';
%     mat_md5 = DataHash(results, Opt);
%     disp(' ')
%     disp('MAT results saved to: ')
%     disp(results)
%     disp(['File MD5 hash: ' mat_md5])
%     
%     % d) Check if physiologically possible, if not, remove
%     %************************
%     % checkind = find(x(:,1) < 0);
%     % x(checkind,:) = [];
%     % tumind(checkind) = [];
%     
%     % ve > 1
%     % checkind = find(x(:,2) > 1);
%     % x(checkind,:) = [];
%     % tumind(checkind) = [];
%     %
%     % % ve < 0
%     % checkind = find(x(:,2) < 0);
%     % x(checkind,:) = [];
%     % tumind(checkind) = [];
%     %
%     % % vp > 1
%     % checkind = find(x(:,3) > 1);
%     % x(checkind,:) = [];
%     % tumind(checkind) = [];
%     %
%     % % vp < 0
%     % checkind = find(x(:,3) < 0);
%     % x(checkind,:) = [];
%     % tumind(checkind) = [];
%     %
%     % % Ktrans > 3
%     % checkind = find(x(:,1) > 3);
%     % x(checkind,:) = [];
%     % tumind(checkind) = [];
%     %
%     % % Ktrans < 0
%     % checkind = find(x(:,1) < 0);
%     % x(checkind,:) = [];
%     % tumind(checkind) = [];
%     %
%     % %% d.2) Check R2 fit
%     %
%     % checkind = find(x(:,4) < r2filter);
%     % x(checkind,:)    = [];
%     % tumind(checkind) = [];
%     
%     
%     % f) Make maps and Save image files
%     %************************
%     fit_voxels = 0; %% DEBUG
%     %[discard, actual] = fileparts(strrep(dynam_name, '\', '/'));
%     
%     if strcmp(cur_dce_model, 'aif')
%         dce_model_name = 'aif';
%         % Write ROI results
%         if number_rois~=0
%             headings = {'ROI path', 'ROI', 'Ktrans', 'Ve','Residual', 'Ktrans 95% low', ...
%                 'Ktrans 95% high', 'Ve 95% low', 'Ve 95% high'};
%             roi_results(:,3) = []; %erase vp column
%             xls_results = [roi_list roi_name mat2cell(roi_results,ones(1,size(roi_results,1)),ones(1,size(roi_results,2)))];
%             xls_results = [headings; xls_results];
%             xls_path = fullfile(PathName, [rootname dce_model_name '_rois.xls']);
%             xlswrite(xls_path,xls_results);
%         end
%         % Write voxel results
%         if fit_voxels
%             
%             KtransROI = zeros(size(currentimg));
%             veROI     = zeros(size(currentimg));
%             residual  = zeros(size(currentimg));
%             ci_95_low_ktrans	= zeros(size(currentimg));
%             ci_95_high_ktrans	= zeros(size(currentimg));
%             ci_95_low_ve		= zeros(size(currentimg));
%             ci_95_high_ve		= zeros(size(currentimg));
%             
%             KtransROI(tumind) = fitting_results(:,1);
%             veROI(tumind)     = fitting_results(:,2);
%             residual(tumind)  = fitting_results(:,4);
%             ci_95_low_ktrans(tumind)	= fitting_results(:,5);
%             ci_95_high_ktrans(tumind)	= fitting_results(:,6);
%             ci_95_low_ve(tumind)		= fitting_results(:,7);
%             ci_95_high_ve(tumind)		= fitting_results(:,8);
%             
%             nii_path{1} = fullfile(PathName, [rootname dce_model_name '_Ktrans.nii']);
%             nii_path{2} = fullfile(PathName, [rootname dce_model_name '_ve.nii']);
%             nii_path{3} = fullfile(PathName, [rootname dce_model_name '_residual.nii']);
%             nii_path{4} = fullfile(PathName, [rootname dce_model_name '_ktrans_ci_low.nii']);
%             nii_path{5} = fullfile(PathName, [rootname dce_model_name '_ktrans_ci_high.nii']);
%             nii_path{6} = fullfile(PathName, [rootname dce_model_name '_ve_ci_low.nii']);
%             nii_path{7} = fullfile(PathName, [rootname dce_model_name '_ve_ci_high.nii']);
%             
%             save_nii(make_nii(KtransROI, res, [1 1 1]), nii_path{1});
%             save_nii(make_nii(veROI, res, [1 1 1]), nii_path{2});
%             save_nii(make_nii(residual, res, [1 1 1]), nii_path{3});
%             save_nii(make_nii(ci_95_low_ktrans, res, [1 1 1]), nii_path{4});
%             save_nii(make_nii(ci_95_high_ktrans, res, [1 1 1]), nii_path{5});
%             save_nii(make_nii(ci_95_low_ve, res, [1 1 1]), nii_path{6});
%             save_nii(make_nii(ci_95_high_ve, res, [1 1 1]), nii_path{7});
%         end
%     elseif strcmp(cur_dce_model, 'aif_vp')
%         dce_model_name = 'aif_vp';
%         % Write ROI results
%         if number_rois~=0
%             headings = {'ROI path', 'ROI', 'Ktrans', 'Ve','Vp','Residual', 'Ktrans 95% low', ...
%                 'Ktrans 95% high', 'Ve 95% low', 'Ve 95% high','Vp 95% low','Vp 95% high'};
%             xls_results = [roi_list roi_name mat2cell(roi_results,ones(1,size(roi_results,1)),ones(1,size(roi_results,2)))];
%             xls_results = [headings; xls_results];
%             xls_path = fullfile(PathName, [rootname dce_model_name '_rois.xls']);
%             xlswrite(xls_path,xls_results);
%         end
%         % Write voxel results
%         if fit_voxels
%             KtransROI = zeros(size(currentimg));
%             veROI     = zeros(size(currentimg));
%             vpROI     = zeros(size(currentimg));
%             residual  = zeros(size(currentimg));
%             ci_95_low_ktrans	= zeros(size(currentimg));
%             ci_95_high_ktrans	= zeros(size(currentimg));
%             ci_95_low_ve		= zeros(size(currentimg));
%             ci_95_high_ve		= zeros(size(currentimg));
%             ci_95_low_vp		= zeros(size(currentimg));
%             ci_95_high_vp		= zeros(size(currentimg));
%             
%             KtransROI(tumind) = fitting_results(:,1);
%             veROI(tumind)     = fitting_results(:,2);
%             vpROI(tumind)     = fitting_results(:,3);
%             residual(tumind)  = fitting_results(:,4);
%             ci_95_low_ktrans(tumind)	= fitting_results(:,5);
%             ci_95_high_ktrans(tumind)	= fitting_results(:,6);
%             ci_95_low_ve(tumind)		= fitting_results(:,7);
%             ci_95_high_ve(tumind)		= fitting_results(:,8);
%             ci_95_low_vp(tumind)		= fitting_results(:,9);
%             ci_95_high_vp(tumind)		= fitting_results(:,10);
%             
%             nii_path{1} = fullfile(PathName, [rootname dce_model_name '_Ktrans.nii']);
%             nii_path{2} = fullfile(PathName, [rootname dce_model_name '_ve.nii']);
%             nii_path{3} = fullfile(PathName, [rootname dce_model_name '_vp.nii']);
%             nii_path{4} = fullfile(PathName, [rootname dce_model_name '_residual.nii']);
%             nii_path{5} = fullfile(PathName, [rootname dce_model_name '_ktrans_ci_low.nii']);
%             nii_path{6} = fullfile(PathName, [rootname dce_model_name '_ktrans_ci_high.nii']);
%             nii_path{7} = fullfile(PathName, [rootname dce_model_name '_ve_ci_low.nii']);
%             nii_path{8} = fullfile(PathName, [rootname dce_model_name '_ve_ci_high.nii']);
%             nii_path{9} = fullfile(PathName, [rootname dce_model_name '_vp_ci_low.nii']);
%             nii_path{10} = fullfile(PathName, [rootname dce_model_name '_vp_ci_high.nii']);
%             
%             save_nii(make_nii(KtransROI, res, [1 1 1]), nii_path{1});
%             save_nii(make_nii(veROI, res, [1 1 1]), nii_path{2});
%             save_nii(make_nii(vpROI, res, [1 1 1]), nii_path{3});
%             save_nii(make_nii(residual, res, [1 1 1]), nii_path{4});
%             save_nii(make_nii(ci_95_low_ktrans, res, [1 1 1]), nii_path{5});
%             save_nii(make_nii(ci_95_high_ktrans, res, [1 1 1]), nii_path{6});
%             save_nii(make_nii(ci_95_low_ve, res, [1 1 1]), nii_path{7});
%             save_nii(make_nii(ci_95_high_ve, res, [1 1 1]), nii_path{8});
%             save_nii(make_nii(ci_95_low_vp, res, [1 1 1]), nii_path{9});
%             save_nii(make_nii(ci_95_high_vp, res, [1 1 1]), nii_path{10});
%         end
%     elseif strcmp(cur_dce_model, 'fxr')
%         dce_model_name = 'fxr';
%         % Write ROI results
%         if number_rois~=0
%             headings = {'ROI path', 'ROI', 'Ktrans', 'Ve','Tau','Residual', 'Ktrans 95% low', ...
%                 'Ktrans 95% high', 'Ve 95% low', 'Ve 95% high','Tau 95% low','Tau 95% high'};
%             xls_results = [roi_list roi_name mat2cell(roi_results,ones(1,size(roi_results,1)),ones(1,size(roi_results,2)))];
%             xls_results = [headings; xls_results];
%             xls_path = fullfile(PathName, [rootname dce_model_name '_rois.xls']);
%             xlswrite(xls_path,xls_results);
%         end
%         % Write voxel results
%         if fit_voxels
%             KtransROI = zeros(size(currentimg));
%             veROI     = zeros(size(currentimg));
%             tauROI     = zeros(size(currentimg));
%             residual  = zeros(size(currentimg));
%             ci_95_low_ktrans	= zeros(size(currentimg));
%             ci_95_high_ktrans	= zeros(size(currentimg));
%             ci_95_low_ve		= zeros(size(currentimg));
%             ci_95_high_ve		= zeros(size(currentimg));
%             ci_95_low_tau		= zeros(size(currentimg));
%             ci_95_high_tau		= zeros(size(currentimg));
%             
%             KtransROI(tumind) = fitting_results(:,1);
%             veROI(tumind)     = fitting_results(:,2);
%             tauROI(tumind)     = fitting_results(:,3);
%             residual(tumind)  = fitting_results(:,4);
%             ci_95_low_ktrans(tumind)	= fitting_results(:,5);
%             ci_95_high_ktrans(tumind)	= fitting_results(:,6);
%             ci_95_low_ve(tumind)		= fitting_results(:,7);
%             ci_95_high_ve(tumind)		= fitting_results(:,8);
%             ci_95_low_tau(tumind)		= fitting_results(:,9);
%             ci_95_high_tau(tumind)		= fitting_results(:,10);
%             
%             nii_path{1} = fullfile(PathName, [rootname dce_model_name '_Ktrans.nii']);
%             nii_path{2} = fullfile(PathName, [rootname dce_model_name '_ve.nii']);
%             nii_path{3} = fullfile(PathName, [rootname dce_model_name '_tau.nii']);
%             nii_path{4} = fullfile(PathName, [rootname dce_model_name '_residual.nii']);
%             nii_path{5} = fullfile(PathName, [rootname dce_model_name '_ktrans_ci_low.nii']);
%             nii_path{6} = fullfile(PathName, [rootname dce_model_name '_ktrans_ci_high.nii']);
%             nii_path{7} = fullfile(PathName, [rootname dce_model_name '_ve_ci_low.nii']);
%             nii_path{8} = fullfile(PathName, [rootname dce_model_name '_ve_ci_high.nii']);
%             nii_path{9} = fullfile(PathName, [rootname dce_model_name '_tau_ci_low.nii']);
%             nii_path{10} = fullfile(PathName, [rootname dce_model_name '_tau_ci_high.nii']);
%             
%             save_nii(make_nii(KtransROI, res, [1 1 1]), nii_path{1});
%             save_nii(make_nii(veROI, res, [1 1 1]), nii_path{2});
%             save_nii(make_nii(tauROI, res, [1 1 1]), nii_path{3});
%             save_nii(make_nii(residual, res, [1 1 1]), nii_path{4});
%             save_nii(make_nii(ci_95_low_ktrans, res, [1 1 1]), nii_path{5});
%             save_nii(make_nii(ci_95_high_ktrans, res, [1 1 1]), nii_path{6});
%             save_nii(make_nii(ci_95_low_ve, res, [1 1 1]), nii_path{7});
%             save_nii(make_nii(ci_95_high_ve, res, [1 1 1]), nii_path{8});
%             save_nii(make_nii(ci_95_low_tau, res, [1 1 1]), nii_path{9});
%             save_nii(make_nii(ci_95_high_tau, res, [1 1 1]), nii_path{10});
%         end
%     elseif strcmp(cur_dce_model, 'fractal')
%     elseif strcmp(cur_dce_model, 'auc')
%         
%         dce_model_name = 'auc';
%         % Write ROI results
%         if number_rois~=0
%             headings = {'ROI path', 'ROI', 'AUC conc', 'AUC sig','NAUC conc', 'NAUC sig'};
%    
%             xls_results = [roi_list roi_name mat2cell(roi_results,ones(1,size(roi_results,1)),ones(1,size(roi_results,2)))];
%             xls_results = [headings; xls_results];
%             xls_path = fullfile(PathName, [rootname dce_model_name '_rois.xls']);
%             xlswrite(xls_path,xls_results);
%         end
%         
%          % Write voxel results
%         if fit_voxels
%             AUCc     = zeros(size(currentimg));
%             AUCs     = zeros(size(currentimg));
%             NAUCc    = zeros(size(currentimg));
%             NAUCs    = zeros(size(currentimg));
%             
%             
%             AUCc(tumind)     = fitting_results(:,1);
%             AUCs(tumind)     = fitting_results(:,2);
%             NAUCc(tumind)    = fitting_results(:,3);
%             NAUCs(tumind)    = fitting_results(:,4);
%            
%             nii_path{1} = fullfile(PathName, [rootname dce_model_name '_AUCc.nii']);
%             nii_path{2} = fullfile(PathName, [rootname dce_model_name '_AUCs.nii']);
%             nii_path{3} = fullfile(PathName, [rootname dce_model_name '_NAUCs.nii']);
%             nii_path{4} = fullfile(PathName, [rootname dce_model_name '_NAUCs.nii']);
%             
%             
%             save_nii(make_nii(AUCs, res, [1 1 1]), nii_path{1});
%             save_nii(make_nii(AUCs, res, [1 1 1]), nii_path{2});
%             save_nii(make_nii(NAUCc, res, [1 1 1]), nii_path{3});
%             save_nii(make_nii(NAUCs, res, [1 1 1]), nii_path{4});
%            
%         end
%         
%         
% 
%     else
%         error('Model not supported');
%     end
%     
%     % Calculate file hashes and log them
%     Opt.Input = 'file';
%     if number_rois~=0
%         xls_md5 = DataHash(xls_path, Opt);
%         disp(' ')
%         disp('ROI results saved to: ')
%         disp(xls_path)
%         disp(['File MD5 hash: ' xls_md5])
%     end
%     if fit_voxels
%         disp(' ')
%         disp('Voxel results saved to: ')
%         for i=1:numel(nii_path)
%             nii_md5 = DataHash(nii_path{i}, Opt);
%             disp(nii_path{i})
%             disp(['File MD5 hash: ' nii_md5])
%         end
%     end
%     
%     disp(['Finished making maps for ' cur_dce_model '...']);
%     
%  end
    
    disp(' ');
    disp('Finished D');
    disp(datestr(now))
    toc
    diary off;
    



