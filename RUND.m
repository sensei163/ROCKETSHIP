function varargout = RUND(varargin)
% RUND MATLAB code for RUND.fig
%      RUND, by itself, creates a new RUND or raises the existing
%      singleton*.
%
%      H = RUND returns the handle to a new RUND or the handle to
%      the existing singleton*.
%
%      RUND('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RUND.M with the given input arguments.
%
%      RUND('Property','Value',...) creates a new RUND or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RUND_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RUND_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RUND

% Last Modified by GUIDE v2.5 03-Feb-2014 17:50:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RUND_OpeningFcn, ...
                   'gui_OutputFcn',  @RUND_OutputFcn, ...
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


% --- Executes just before RUND is made visible.
function RUND_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RUND (see VARARGIN)

% Choose default command line output for RUND
handles.output = hObject;

% Setup the name of the B results file
BRESULTS = varargin{1};

set(handles.results_b_path, 'String', BRESULTS{1});

handles.batch = 0; % Batch central 

% Update handles structure
guidata(hObject, handles);

uirestore;

uirestore(handles.fractal);
uirestore(handles.auc);
uirestore(handles.fxr);
uirestore(handles.aif);
uirestore(handles.aif_vp);
uirestore(handles.none);
uirestore(handles.moving);
uirestore(handles.rlowess);

% check preference file withe verbose for invalid selections, use this to 
% guard against typos
parse_preference_file('dce_preferences.txt',1,...
    {'aif_lower_limits' 'aif_upper_limits' 'aif_initial_values' ...
    'aif_TolFun' 'aif_TolX' 'aif_MaxIter' 'aif_MaxFunEvals' 'aif_Robust'...
    'voxel_lower_limit_ktrans' 'voxel_upper_limit_ktrans' 'voxel_initial_value_ktrans' ...
    'voxel_lower_limit_ve' 'voxel_upper_limit_ve' 'voxel_initial_value_ve' ...
    'voxel_lower_limit_vp' 'voxel_upper_limit_vp' 'voxel_initial_value_vp' ...
    'voxel_TolFun' 'voxel_TolX' 'voxel_MaxIter' 'voxel_MaxFunEvals' ...
    'voxel_lower_limit_tau' 'voxel_upper_limit_tau' 'voxel_initial_value_tau' ...
    'voxel_Robust' 'fxr_fw'});
% Create structure to hold roi list
handles.roi_list = {};
handles.saved_results = '';

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RUND wait for user response (see UIRESUME)
 uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RUND_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} =handles.saved_results;
varargout{2} = handles.batch;
delete(handles.figure1);



function results_b_path_Callback(hObject, eventdata, handles)
% hObject    handle to results_b_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of results_b_path as text
%        str2double(get(hObject,'String')) returns contents of results_b_path as a double


% --- Executes during object creation, after setting all properties.
function results_b_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to results_b_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_results_b.
function browse_results_b_Callback(hObject, eventdata, handles)
% hObject    handle to browse_results_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);

[filename, pathname, filterindex] = uigetfile( ...
    {  '*.mat','Matlab Worksapce Files (*.mat)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Choose Results from Section B'); %#ok<NASGU>
if isequal(filename,0)
    %disp('User selected Cancel')
else
    %disp(['User selected ', fullfile(pathname, filename)])
    
    % Combine path and filename together
    fullpath = strcat(pathname,filename);

    set(handles.results_b_path,'String',fullpath);
end
guidata(hObject, handles);



function time_smoothing_window_Callback(hObject, eventdata, handles)
% hObject    handle to time_smoothing_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_smoothing_window as text
%        str2double(get(hObject,'String')) returns contents of time_smoothing_window as a double
uiremember;

% --- Executes during object creation, after setting all properties.
function time_smoothing_window_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_smoothing_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xy_smooth_size_Callback(hObject, eventdata, handles)
% hObject    handle to xy_smooth_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xy_smooth_size as text
%        str2double(get(hObject,'String')) returns contents of xy_smooth_size as a double
uiremember;

% --- Executes during object creation, after setting all properties.
function xy_smooth_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xy_smooth_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in roi_box.
function roi_box_Callback(hObject, eventdata, handles)
% hObject    handle to roi_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns roi_box contents as cell array
%        contents{get(hObject,'Value')} returns selected item from roi_box


% --- Executes during object creation, after setting all properties.
function roi_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roi_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in add_roi.
function add_roi_Callback(hObject, eventdata, handles)
% hObject    handle to add_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);

[filename, pathname, filterindex] = uigetfile( ...
    {  '*.nii','Nifti Files (*.nii)'; ...
    '*.hdr;*.img','Analyze Files (*.hdr, *.img)';...
    '*.*',  'All Files (*.*)'}, ...
    'Pick a file', ...
    'MultiSelect', 'on'); %#ok<NASGU>
if isequal(filename,0)
    %disp('User selected Cancel')
else
    %disp(['User selected ', fullfile(pathname, filename)])
    list = get(handles.roi_box,'String');
    
    % Combine path and filename together
    fullpath = strcat(pathname,filename);
    
    % Stupid matlab uses a different datastructure if only one file
    % is selected, handle special case
    if ischar(list)
        list = {list};
    end
    if ischar(filename)
        filename = {filename};
    end
    if ischar(fullpath)
        fullpath = {fullpath};
    end

    filename = filename';
    fullpath = fullpath';
        
    % Add selected files to listbox
    if strcmp(list,'No Files')
        list = filename;
        handles.roi_list = fullpath;
    else
        list = [list;  filename];
        handles.roi_list = [handles.roi_list; fullpath];
    end 
    
    set(handles.roi_box,'String',list, 'Value',1)
end
guidata(hObject, handles);

% --- Executes on button press in remove_roi.
function remove_roi_Callback(hObject, eventdata, handles)
% hObject    handle to remove_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index_selected = get(handles.roi_box,'Value');
list = get(handles.roi_box,'String');
for n=size(index_selected,2):-1:1
    % Remove from end of list first so resizing does not 
    % change subsequent index numbers
    %disp(['User removed ', list{index_selected(n)}]);
    list(index_selected(n)) = [];
    handles.roi_list(index_selected(n)) = [];
end

set(handles.roi_box,'String',list, 'Value',1)
guidata(hObject, handles);

% --- Executes on button press in fit_voxels.
function fit_voxels_Callback(hObject, eventdata, handles)
% hObject    handle to fit_voxels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fit_voxels
uiremember;


function number_cpus_Callback(hObject, eventdata, handles)
% hObject    handle to number_cpus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_cpus as text
%        str2double(get(hObject,'String')) returns contents of number_cpus as a double


% --- Executes during object creation, after setting all properties.
function number_cpus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_cpus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in neuroecon.
function neuroecon_Callback(hObject, eventdata, handles)
% hObject    handle to neuroecon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of neuroecon


% --- Executes on button press in done.
function done_Callback(hObject, eventdata, handles)
% hObject    handle to done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
results_b_path = get(handles.results_b_path,'String');

dce_model.aif    = get(handles.aif, 'Value');
dce_model.aif_vp = get(handles.aif_vp, 'Value');
dce_model.fxr    = get(handles.fxr, 'Value');

dce_model.fractal= get(handles.fractal, 'Value');
dce_model.auc    = get(handles.auc, 'Value');


time_smoothing = get(get(handles.time_smoothing,'SelectedObject'),'Tag');
time_smoothing_window = str2num(get(handles.time_smoothing_window, 'String')); %#ok<ST2NM>
xy_smooth_size = str2num(get(handles.xy_smooth_size, 'String')); %#ok<ST2NM>
number_cpus = str2num(get(handles.number_cpus, 'String')); %#ok<ST2NM>
neuroecon = get(handles.neuroecon, 'Value'); 

roi_list = handles.roi_list;
fit_voxels = get(handles.fit_voxels,'Value');
batch  = handles.batch;
outputft = get(handles.outputft, 'Value');

saved_results = D_fit_voxels_func(results_b_path,dce_model,time_smoothing,time_smoothing_window,xy_smooth_size,number_cpus,roi_list,fit_voxels,neuroecon, batch, outputft);

handles.saved_results = saved_results;
handles.batch = 0;
guidata(hObject, handles);

uiresume(handles.figure1);



% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
% hObject    handle to cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargout{1} = ''; %handles.output;
% delete(handles.figure1);
uiresume(handles.figure1);

% --- Executes on button press in aif.
function aif_Callback(hObject, eventdata, handles)
% hObject    handle to aif (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of aif

%set(handles.aif_vp, 'Value', 0);
uiremember;


% --- Executes on button press in aif_vp.
function aif_vp_Callback(hObject, eventdata, handles)
% hObject    handle to aif_vp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%set(handles.aif, 'Value', 0);
uiremember;
% Hint: get(hObject,'Value') returns toggle state of aif_vp


% --- Executes on button press in fxr.
function fxr_Callback(hObject, eventdata, handles)
% hObject    handle to fxr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiremember;
% Hint: get(hObject,'Value') returns toggle state of fxr


% --- Executes on button press in sauc.
function sauc_Callback(hObject, eventdata, handles)
% hObject    handle to sauc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiremember;
% Hint: get(hObject,'Value') returns toggle state of sauc


% --- Executes on button press in fractal.
function fractal_Callback(hObject, eventdata, handles)
% hObject    handle to fractal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiremember;
% Hint: get(hObject,'Value') returns toggle state of fractal


% --- Executes on button press in auc.
function auc_Callback(hObject, eventdata, handles)
% hObject    handle to auc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiremember;
% Hint: get(hObject,'Value') returns toggle state of auc


% --- Executes on button press in auc_rr.
function auc_rr_Callback(hObject, eventdata, handles)
% hObject    handle to auc_rr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiremember;
% Hint: get(hObject,'Value') returns toggle state of auc_rr


% --- Executes during object deletion, before destroying properties.
function aif_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to aif (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over aif.
function aif_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to aif (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on aif and none of its controls.
function aif_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to aif (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in none.
function none_Callback(hObject, eventdata, handles)
% hObject    handle to none (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of none


% --- Executes when selected object is changed in time_smoothing.
function time_smoothing_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in time_smoothing 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
uiremember;


% --- Executes during object creation, after setting all properties.
function aif_CreateFcn(hObject, eventdata, handles)
% hObject    handle to aif (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in prep_batch.
function prep_batch_Callback(hObject, eventdata, handles)
% hObject    handle to prep_batch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

results_b_path = get(handles.results_b_path,'String');

dce_model.aif    = get(handles.aif, 'Value');
dce_model.aif_vp = get(handles.aif_vp, 'Value');
dce_model.fxr    = get(handles.fxr, 'Value');
%dce_model.sauc   = get(handles.sauc, 'Value');
dce_model.fractal= get(handles.fractal, 'Value');
dce_model.auc    = get(handles.auc, 'Value');
%dce_model.auc_rr = get(handles.auc_rr, 'Value');

time_smoothing = get(get(handles.time_smoothing,'SelectedObject'),'Tag');
time_smoothing_window = str2num(get(handles.time_smoothing_window, 'String')); %#ok<ST2NM>
xy_smooth_size = str2num(get(handles.xy_smooth_size, 'String')); %#ok<ST2NM>
number_cpus = str2num(get(handles.number_cpus, 'String')); %#ok<ST2NM>
neuroecon = get(handles.neuroecon, 'Value'); 

roi_list = handles.roi_list;
fit_voxels = get(handles.fit_voxels,'Value');
batch  = 1;
outputft = get(handles.outputft, 'Value');
[saved_results, batch] = D_fit_voxels_func(results_b_path,dce_model,time_smoothing,time_smoothing_window,xy_smooth_size,number_cpus,roi_list,fit_voxels,neuroecon, batch, outputft);

handles.saved_results = saved_results;
handles.batch         = 1;
uiremember;
guidata(hObject, handles);
uiresume(handles.figure1);


% --- Executes on button press in FXL_rr.
function FXL_rr_Callback(hObject, eventdata, handles)
% hObject    handle to FXL_rr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FXL_rr


% --- Executes on button press in FXR_rr.
function FXR_rr_Callback(hObject, eventdata, handles)
% hObject    handle to FXR_rr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FXR_rr


% --- Executes on selection change in outputft.
function outputft_Callback(hObject, eventdata, handles)
% hObject    handle to outputft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns outputft contents as cell array
%        contents{get(hObject,'Value')} returns selected item from outputft


% --- Executes during object creation, after setting all properties.
function outputft_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outputft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
% delete(hObject);
uiresume(hObject);
