function varargout = fitting_analysis(varargin)
% FITTING_ANALYSIS MATLAB code for fitting_analysis.fig
%      FITTING_ANALYSIS, by itself, creates a new FITTING_ANALYSIS or raises the existing
%      singleton*.
%
%      H = FITTING_ANALYSIS returns the handle to a new FITTING_ANALYSIS or the handle to
%      the existing singleton*.
%
%      FITTING_ANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FITTING_ANALYSIS.M with the given input arguments.
%
%      FITTING_ANALYSIS('Property','Value',...) creates a new FITTING_ANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fitting_analysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fitting_analysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fitting_analysis

% Last Modified by GUIDE v2.5 28-Jan-2014 10:34:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fitting_analysis_OpeningFcn, ...
                   'gui_OutputFcn',  @fitting_analysis_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1}) && ~strcmp(varargin{1},'results_path')
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before fitting_analysis is made visible.
function fitting_analysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fitting_analysis (see VARARGIN)

% Choose default command line output for fitting_analysis
handles.output = hObject;
handles.roi_data_ready = 0;
handles.voxel_data_ready = 0;

% Get function inputs
if nargin>1 && numel(varargin)>1 && strcmp(varargin{1},'results_path')
    set(handles.results_cfit_path,'String',varargin{2});
    % Update cfit structures
    handles = cfit_path_changed(handles);
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fitting_analysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = fitting_analysis_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


function results_cfit_path_Callback(hObject, eventdata, handles)
handles = cfit_path_changed(handles);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function results_cfit_path_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_cfit_results.
function browse_cfit_results_Callback(hObject, eventdata, handles)
guidata(hObject, handles);

[filename, pathname, filterindex] = uigetfile( ...
    {  '*.mat','Matlab Worksapce Files (*.mat)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Choose Curve Fitting Results'); %#ok<NASGU>
if isequal(filename,0)
    %disp('User selected Cancel')
else
    %disp(['User selected ', fullfile(pathname, filename)])
    
    % Combine path and filename together
    fullpath = strcat(pathname,filename);

    set(handles.results_cfit_path,'String',fullpath);
end

% Update cfit structures
handles = cfit_path_changed(handles);

% Update handles structure
guidata(hObject, handles);


function background_image_path_Callback(hObject, eventdata, handles)
handles = load_check_data(handles);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function background_image_path_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_image.
function browse_image_Callback(hObject, eventdata, handles)
guidata(hObject, handles);

[filename, pathname, filterindex] = uigetfile( ...
    {  '*.nii','Nifti Files (*.nii)'; ...
    '*2dseq','Bruker Files (2dseq)'; ...
    '*.hdr;*.img','Analyze Files (*.hdr, *.img)';...
    '*.*',  'All Files (*.*)'}, ...
    'Choose Background Image file'); %#ok<NASGU>
if isequal(filename,0)
    %disp('User selected Cancel')
else
    %disp(['User selected ', fullfile(pathname, filename)])
    
    % Combine path and filename together
    fullpath = strcat(pathname,filename);

    set(handles.background_image_path,'String',fullpath);
end

handles = load_check_data(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_run.
function button_run_Callback(hObject, eventdata, handles)
if handles.voxel_data_ready
    results_cfit_path = get(handles.results_cfit_path,'String');
    background_image_path = get(handles.background_image_path,'String');
    show_original = get(handles.show_original,'Value');
    show_ci = get(handles.show_ci,'Value');

    compare_fits(results_cfit_path,background_image_path,show_original,show_ci);
else
    handles = load_check_data(handles,'voxel');
    % Update handles structure
    guidata(hObject, handles);
end



% --- Executes on button press in show_original.
function show_original_Callback(hObject, eventdata, handles)


% --- Executes on button press in show_ci.
function show_ci_Callback(hObject, eventdata, handles)


% --- Executes on selection change in roi_listbox.
function roi_listbox_Callback(hObject, eventdata, handles)
if handles.roi_data_ready
    contents = cellstr(get(hObject,'String')); 
    selected_name = contents{get(hObject,'Value')};
    selected_roi_temp = strfind(handles.fit_data.roi_name,selected_name);
    selected_roi = find(not(cellfun('isempty', selected_roi_temp)));

    plot_data.Ct			= handles.xdata.roi_series(:,selected_roi);
    plot_data.Ct_original	= handles.xdata.roi_series_original(:,selected_roi);
    plot_data.Cp			= handles.xdata.Cp;
    plot_data.timer			= handles.xdata.timer;
    plot_data.x_units = 'Time (minutes)';
    plot_data.y_units =  'Concentration (mmol)';
    plot_data.fit_parameters= handles.fit_data.roi_results(selected_roi,:);
    plot_data.model_name		= handles.fit_data.model_name;
    plot_data.show_original = get(handles.show_original,'Value');
    plot_data.show_ci		= get(handles.show_ci,'Value');
    plot_data.title = ['ROI "' selected_name '"'];

    if strcmp(handles.fit_data.model_name,'fxr')
        plot_data.R1o = handles.xdata.roi_r1(selected_roi);
        plot_data.R1i = handles.xdata.roi_r1(selected_roi);
        plot_data.r1 = handles.xdata.relaxivity;
        plot_data.fw = 0.8;
    end

    figure(2);
    plot_dce_curve(plot_data);
else
    handles = load_check_data(handles,'roi');
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function roi_listbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Should be run when path is changed to read in new file and update variables.
function return_handles = cfit_path_changed(handles)
handles = load_check_data(handles);

if handles.roi_data_ready
    set(handles.roi_listbox,'String',handles.fit_data.roi_name, 'Value',1)
else
    set(handles.roi_listbox,'String','', 'Value',0)
end

% Update handles structure
return_handles = handles;



function return_handles = load_check_data(handles,verbose)
if nargin<2
    verbose = 'default';
end
results_cfit_path = get(handles.results_cfit_path,'String');
background_image_path = get(handles.background_image_path,'String');
lower_model_path = get(handles.lower_model_path,'String');
handles.roi_data_ready = 0;
handles.voxel_data_ready = 0;
handles.ftest_ready = 0;
voxel_message = 'No results file selected';
roi_message = 'No results file selected';
ftest_message = 'No primary results file selected';

% Load data and check if ready for voxel and ROI analysis
try
    load(results_cfit_path);
    handles.xdata = xdata{1};
    handles.fit_data = fit_data;
    
    [~, temp_name, temp_ext] = fileparts(results_cfit_path);
    ready_message = ['File loaded: ' temp_name temp_ext];   
    
    if handles.fit_data.number_rois~=0
        handles.roi_data_ready = 1;
    else
        roi_message = 'No ROI results in file';
    end
    if handles.fit_data.fit_voxels 
        if exist(background_image_path,'file')
            handles.voxel_data_ready = 1;
        else
            voxel_message = 'Selected background image does not exist';
        end
    else
        voxel_message = 'No voxel results in file';
    end
        
    information_string = {['Fit Model: ' handles.fit_data.model_name],...
        ['Fitted ROIs: ' num2str(handles.fit_data.number_rois)],...
        ['Fitted Voxels: ' num2str(handles.xdata.numvoxels*handles.fit_data.fit_voxels)]};
    
catch err
    information_string = {};
    handles.roi_data_ready = 0;
    handles.voxel_data_ready = 0;
    if ~exist(results_cfit_path,'file')
        ready_message = 'results file not found';
    elseif ~exist('xdata','var') || ~exist('fit_data','var')
        ready_message = 'results file does not contain fit data';
    elseif ~isfield(handles.fit_data,'fit_voxels') ||...
            ~isfield(handles.fit_data,'number_rois') ||...
            ~isfield(handles.fit_data,'model_name')
        ready_message = 'fit data is incomplete, rerun fit';
    else
        rethrow(err);
    end
    voxel_message = ready_message;
    roi_message = ready_message;
    ftest_message = ready_message;
end

% Now setup for ftest
try
    lower_model = load(lower_model_path); 
    handles.lower_model_xdata = lower_model.xdata{1};
    handles.lower_model_fit_data = lower_model.fit_data;
    
    handles.ftest_ready = 1;
catch err
    handles.ftest_ready = 0;
    if ~exist(lower_model_path,'file')
        ftest_message = 'comparison model file not selected';
        if ~isempty(lower_model_path)
            ftest_message = 'comparison model file not found';
        end
    elseif ~isfield(lower_model,'xdata') || ~isfield(lower_model,'fit_data')
        ftest_message = 'comparison model file does not contain fit data';
    else
        rethrow(err);
    end
    if ~isempty(lower_model_path)
        ready_message = ftest_message;
    end
end

% Display info messages
set(handles.cfit_information,'String',information_string);

if strcmp(verbose,'default')
    if handles.voxel_data_ready || handles.roi_data_ready
        update_status(handles,ready_message, 'black');
    else
        update_status(handles,ready_message, 'red');
    end
elseif strcmp(verbose,'voxel')
    if handles.voxel_data_ready
        update_status(handles,voxel_message, 'black');
    else
        update_status(handles,voxel_message, 'red');
    end
elseif strcmp(verbose,'roi')
    if handles.roi_data_ready
        update_status(handles,roi_message, 'black');
    else
        update_status(handles,roi_message, 'red');
    end
elseif strcmp(verbose,'ftest')
    if handles.ftest_ready
        update_status(handles,ftest_message, 'black');
    else
        update_status(handles,ftest_message, 'red');
    end
end

% Update handles structure
return_handles = handles;



function lower_model_path_Callback(hObject, eventdata, handles)
% hObject    handle to lower_model_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lower_model_path as text
%        str2double(get(hObject,'String')) returns contents of lower_model_path as a double


% --- Executes during object creation, after setting all properties.
function lower_model_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lower_model_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_lower_model.
function browse_lower_model_Callback(hObject, eventdata, handles)
guidata(hObject, handles);

[filename, pathname, filterindex] = uigetfile( ...
    {  '*.mat','Matlab Worksapce Files (*.mat)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Choose Fitting Results'); %#ok<NASGU>
if isequal(filename,0)
    %disp('User selected Cancel')
else
    %disp(['User selected ', fullfile(pathname, filename)])
    
    % Combine path and filename together
    fullpath = strcat(pathname,filename);

    set(handles.lower_model_path,'String',fullpath);
end

% Update structures
handles = load_check_data(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in button_akaike.
function button_akaike_Callback(hObject, eventdata, handles)
run_comparison(handles,2)

% --- Executes on button press in button_ftest.
function button_ftest_Callback(hObject, eventdata, handles)
run_comparison(handles,1)

function run_comparison(handles,test_index)
if handles.ftest_ready
    compare_voxels = (handles.lower_model_fit_data.fit_voxels && handles.fit_data.fit_voxels);
    compare_rois = (handles.lower_model_fit_data.number_rois>0 && handles.fit_data.number_rois>0);
    information_string = get(handles.cfit_information,'String');
    information_string = information_string(1:3);
    information_string(end+1) = {['Compare Model: ' handles.lower_model_fit_data.model_name]};
    % Create custum strings
    if test_index==1
        stat_name = 'p value';
        path_suffix = '_ftest';
        test_name = 'f-test';
    elseif test_index==2
        stat_name = 'relative likelihood';
        path_suffix = '_aic';
        test_name = 'Akaike information criteria';
    end
        
    if compare_voxels 
        disp(['Starting ' test_name ' on voxels']);
        [sse_lower,fp_lower,sse_higher,fp_higher,n]=...
            get_sse_and_fp(handles,1);
        
        % Sanity check
        if fp_lower>=fp_higher
            update_status(handles,'lower model must be a lower number of free parameters','red');
            disp('stopping');
            return;
        end
        
        % Run Test
        number_voxels = numel(sse_higher);
        stat_voxels = 2.*ones(number_voxels,1);
        for i=1:number_voxels
            if test_index==1
                [ p, Fstat, df1, df2 ] = ftest(n,fp_lower,...
                    fp_higher,sse_lower(i),sse_higher(i));
                stat_voxels(i) = p;
            elseif test_index==2
                aic_lower = n*log(sse_lower(i)/n)+2*fp_lower;
                aic_higher = n*log(sse_higher(i)/n)+2*fp_higher;
                relative_likelihood = exp((aic_lower-aic_higher)/2);
                if relative_likelihood>1
                    relative_likelihood = -1+1/relative_likelihood;
                else
                    relative_likelihood = 1-relative_likelihood;
                end
                % positive, relative likelihood lower model minimizes information loss
                % negative, relative likelihood higher model minimizes information loss
                % near +1 lower model better, near -1 higher model better near
                % zero poor inference
                stat_voxels(i) = relative_likelihood;
            end
        end
        
        mean_stat = mean(stat_voxels);  
        disp(['Average voxel ' stat_name ' = ' num2str(mean_stat)]);
        information_string(end+1) = {['Average voxel ' stat_name ' = ' num2str(mean_stat)]};
        
        % Save results
        results_cfit_path = get(handles.results_cfit_path,'String');
        [base_path, ~, ~] = fileparts(results_cfit_path);
        [~, base_name, ~] = fileparts(handles.fit_data.dynam_name);
        save_path = fullfile(base_path,[base_name '_' ...
            handles.fit_data.model_name '_' ...
            handles.lower_model_fit_data.model_name path_suffix '.nii']);
    
        stat_matrix     = zeros([256 256]);
        stat_matrix(handles.fit_data.tumind) = stat_voxels;
        save_nii(make_nii(stat_matrix, [1 1 1], [1 1 1]), save_path);
        disp(['Completed ' test_name ' on voxels']);
    end
    if compare_rois
        disp(['Starting ' test_name ' on ROIs']);
        [sse_lower,fp_lower,sse_higher,fp_higher,n]=...
            get_sse_and_fp(handles,2);
        
        % Sanity check
        if fp_lower>=fp_higher
            update_status(handles,'lower model must be a lower number of free parameters','red');
            disp('stopping');
            return;
        end
        
        % Run Test
        number_rois = numel(sse_higher);
        stat_rois = 2.*ones(number_rois,1);
        f_rois = 2.*ones(number_rois,1);
        for i=1:number_rois
            if test_index==1
                [ p, Fstat, df1, df2 ] = ftest(n,fp_lower,...
                    fp_higher,sse_lower(i),sse_higher(i));
                stat_rois(i) = p;
                f_rois(i) = Fstat;
            elseif test_index==2
                aic_lower = n*log(sse_lower(i)/n)+2*fp_lower;
                aic_higher = n*log(sse_higher(i)/n)+2*fp_higher;
                relative_likelihood = exp((aic_lower-aic_higher)/2);
                if relative_likelihood>1
                    relative_likelihood = -1+1/relative_likelihood;
                else
                    relative_likelihood = 1-relative_likelihood;
                end
                % positive, relative likelihood lower model minimizes information loss
                % negative, relative likelihood higher model minimizes information loss
                % near +1 lower model better, near -1 higher model better near
                % zero poor inference
                stat_rois(i) = relative_likelihood;
            end
        end
        mean_stat = mean(stat_rois);  
        disp(['Average ROI ' stat_name ' = ' num2str(mean_stat)]);
        information_string(end+1) = {['Average ROI ' stat_name ' = ' num2str(mean_stat)]};
        
        % Save results
        results_cfit_path = get(handles.results_cfit_path,'String');
        [base_path, ~, ~] = fileparts(results_cfit_path);
        [~, base_name, ~] = fileparts(handles.fit_data.dynam_name);
        save_path = fullfile(base_path,[base_name '_' ...
            handles.fit_data.model_name '_' ...
            handles.lower_model_fit_data.model_name path_suffix '.xls']);
    
        headings = {'ROI', stat_name, ['Residual ' handles.fit_data.model_name],...
            ['Residual ' handles.lower_model_fit_data.model_name]};
        xls_results = [handles.fit_data.roi_name num2cell(stat_rois) num2cell(sse_higher) num2cell(sse_lower)];
        xls_results = [headings; xls_results];
        
        xlswrite(save_path,xls_results);
               
        disp(['Completed ' test_name ' on ROIs']);
    end
    
    if ~compare_rois && ~compare_voxels
        update_status(handles,'cannot compare, same regions not fitted','red');
    else
        information_string(end+1) = {['lower ' stat_name ' indicates higher order model is better fit']};
        set(handles.cfit_information,'String',information_string);
        disp(['lower ' stat_name ' indicates higher order model is better fit']);
    end
else
    handles = load_check_data(handles,'ftest');
    % Update handles structure
    guidata(hObject, handles);
end
    
function update_status(handles,status_string,color)
set(handles.ready_display,'String',status_string);
set(handles.ready_display, 'ForegroundColor', color);

function [sse1, fp1, sse2, fp2, n]=get_sse_and_fp(handles,region)
if region==1
    % Voxel results
    sse1 = handles.lower_model_fit_data.fitting_results(:,4);
    sse2 = handles.fit_data.fitting_results(:,4);
elseif region==2
    % ROI results
    sse1 = handles.lower_model_fit_data.roi_results(:,4);
    sse2 = handles.fit_data.roi_results(:,4);
end

if strcmp(handles.lower_model_fit_data.model_name,'aif')
    fp1 = 2;
elseif strcmp(handles.lower_model_fit_data.model_name,'aif_vp')
    fp1 = 3;
elseif strcmp(handles.lower_model_fit_data.model_name,'fxr')
    fp1 = 3;
else
    update_status(handles,'selected model not implemented','red');
    return
end

if strcmp(handles.fit_data.model_name,'aif')
    fp2 = 2;
elseif strcmp(handles.fit_data.model_name,'aif_vp')
    fp2 = 3;
elseif strcmp(handles.fit_data.model_name,'fxr')
    fp2 = 3;
else
    update_status(handles,'selected model not implemented','red');
    return
end

n_lower = numel(handles.lower_model_xdata.timer);
n_higher = numel(handles.xdata.timer);

% Sanity check
if n_lower~=n_higher
    update_status(handles,'fits have different number of points','red');
    return;
end

n=n_higher;
