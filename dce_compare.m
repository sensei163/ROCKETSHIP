function varargout = dce_compare(varargin)
% DCE_COMPARE MATLAB code for dce_compare.fig
%      DCE_COMPARE, by itself, creates a new DCE_COMPARE or raises the existing
%      singleton*.
%
%      H = DCE_COMPARE returns the handle to a new DCE_COMPARE or the handle to
%      the existing singleton*.
%
%      DCE_COMPARE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DCE_COMPARE.M with the given input arguments.
%
%      DCE_COMPARE('Property','Value',...) creates a new DCE_COMPARE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dce_compare_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dce_compare_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dce_compare

% Last Modified by GUIDE v2.5 22-Jan-2014 18:42:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dce_compare_OpeningFcn, ...
                   'gui_OutputFcn',  @dce_compare_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before dce_compare is made visible.
function dce_compare_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dce_compare (see VARARGIN)

% Choose default command line output for dce_compare
handles.output = hObject;
handles.roi_data_ready = 0;
handles.voxel_data_ready = 0;


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dce_compare wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dce_compare_OutputFcn(hObject, eventdata, handles) 
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



function return_handles = load_check_data(handles)
results_cfit_path = get(handles.results_cfit_path,'String');
background_image_path = get(handles.background_image_path,'String');
lower_model_path = get(handles.lower_model_path,'String');
handles.roi_data_ready = 0;
handles.voxel_data_ready = 0;
handles.ftest_ready = 0;

% Load data and check if ready for voxel and ROI analysis
try
    load(results_cfit_path);
    handles.xdata = xdata{1};
    handles.fit_data = fit_data;
    
    [~, temp_name, temp_ext] = fileparts(results_cfit_path);
    ready_message = ['File loaded: ' temp_name temp_ext];   
    
    if handles.fit_data.number_rois~=0
        handles.roi_data_ready = 1;
    end
    if handles.fit_data.fit_voxels 
        if exist(background_image_path,'file')
            handles.voxel_data_ready = 1;
        else
            ready_message = 'Selected background image does not exist';
        end
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
        if ~isempty(lower_model_path)
            ready_message = 'lower model file not found';
        end
    elseif ~exist('lower_model.xdata','var') || ~exist('lower_model.fit_data','var')
        if ~isempty(lower_model_path)
            ready_message = 'lower model file does not contain fit data';
        end
    else
        rethrow(err);
    end
end

% Display info messages
set(handles.cfit_information,'String',information_string);

if handles.voxel_data_ready || handles.roi_data_ready
    update_status(handles,ready_message, 'black');
else
    update_status(handles,ready_message, 'red');
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



% --- Executes on button press in button_ftest.
function button_ftest_Callback(hObject, eventdata, handles)
if handles.ftest_ready
    if handles.lower_model_fit_data.fit_voxels && handles.fit_data.fit_voxels 
        [sse_lower,fp_lower,sse_higher,fp_higher,n]=...
            get_sse_and_fp(handles);
        
        % Sanity check
        if fp_lower>=fp_higher
            update_status(handles,'lower model must be a lower number of free parameters','red');
            return;
        end
        
        % Run Test
        number_voxels = numel(sse_higher);
        p_voxels = 2.*ones(number_voxels,1);
        for i=1:number_voxels
            [ p, Fstat, df1, df2 ] = ftest(n,fp_lower,...
                fp_higher,sse_lower(i),sse_higher(i));
            p_voxels(i) = p;
        end
        mean_p = mean(p_voxels);
        information_string = get(handles.cfit_information,'String');
        information_string(end+1) = {['Average p value = ' num2str(mean_p)]};
        information_string(end+1) = {'lower p value indicates higher order model is better fit'};
        set(handles.cfit_information,'String',information_string);
        
        % Save results
        results_cfit_path = get(handles.results_cfit_path,'String');
        [base_path, ~, ~] = fileparts(results_cfit_path);
        [~, base_name, ~] = fileparts(handles.fit_data.dynam_name);
        save_path = fullfile(base_path,[base_name '_' ...
            handles.fit_data.model_name '_' ...
            handles.lower_model_fit_data.model_name '_ftest.nii']);
    
        p_matrix     = zeros([256 256]);
        p_matrix(handles.fit_data.tumind) = p_voxels;
        save_nii(make_nii(p_matrix, [1 1 1], [1 1 1]), save_path);
    end
end
    
    
function update_status(handles,status_string,color)
set(handles.ready_display,'String',status_string);
set(handles.ready_display, 'ForegroundColor', color);

function [sse1, fp1, sse2, fp2, n]=get_sse_and_fp(handles)
if strcmp(handles.lower_model_fit_data.model_name,'aif')
    sse1 = handles.lower_model_fit_data.fitting_results(:,4);
    fp1 = 2;
elseif strcmp(handles.lower_model_fit_data.model_name,'aif_vp')
    sse1 = handles.lower_model_fit_data.fitting_results(:,4);
    fp1 = 3;
elseif strcmp(handles.lower_model_fit_data.model_name,'fxr')
    sse1 = handles.lower_model_fit_data.fitting_results(:,4);
    fp1 = 3;
else
    update_status(handles,'selected model not implemented','red');
    return
end

if strcmp(handles.fit_data.model_name,'aif')
    sse2 = handles.fit_data.fitting_results(:,4);
    fp2 = 2;
elseif strcmp(handles.fit_data.model_name,'aif_vp')
    sse2 = handles.fit_data.fitting_results(:,4);
    fp2 = 3;
elseif strcmp(handles.fit_data.model_name,'fxr')
    sse2 = handles.fit_data.fitting_results(:,4);
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


% --- Executes on button press in button_akaike.
function button_akaike_Callback(hObject, eventdata, handles)
if handles.ftest_ready
    if handles.lower_model_fit_data.fit_voxels && handles.fit_data.fit_voxels 
        [sse_lower,fp_lower,sse_higher,fp_higher,n]=...
            get_sse_and_fp(handles);
        
        % Run Test
        number_voxels = numel(sse_higher);
        aic_voxels = 2.*ones(number_voxels,1);
        for i=1:number_voxels
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
            aic_voxels(i) = relative_likelihood;
        end
        mean_aic = mean(aic_voxels);
        information_string = get(handles.cfit_information,'String');
        information_string(end+1) = {['Average relative likelihood value = ' num2str(mean_aic)]};
        set(handles.cfit_information,'String',information_string);
        
        % Save results
        results_cfit_path = get(handles.results_cfit_path,'String');
        [base_path, ~, ~] = fileparts(results_cfit_path);
        [~, base_name, ~] = fileparts(handles.fit_data.dynam_name);
        save_path = fullfile(base_path,[base_name '_' ...
            handles.fit_data.model_name '_' ...
            handles.lower_model_fit_data.model_name '_aic.nii']);
    
        aic_matrix     = zeros([256 256]);
        aic_matrix(handles.fit_data.tumind) = aic_voxels;
        save_nii(make_nii(aic_matrix, [1 1 1], [1 1 1]), save_path);
    end
end
