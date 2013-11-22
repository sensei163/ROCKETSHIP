function varargout = dce(varargin)
% DCE MATLAB code for dce.fig
%      DCE, by itself, creates a new DCE or raises the existing
%      singleton*.
%
%      H = DCE returns the handle to a new DCE or the handle to
%      the existing singleton*.
%
%      DCE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DCE.M with the given input arguments.
%
%      DCE('Property','Value',...) creates a new DCE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dce_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dce_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dce

% Last Modified by GUIDE v2.5 21-Nov-2013 17:14:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dce_OpeningFcn, ...
                   'gui_OutputFcn',  @dce_OutputFcn, ...
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


% --- Executes just before dce is made visible.
function dce_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dce (see VARARGIN)

% Choose default command line output for dce
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dce wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dce_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in browse_dce.
function browse_dce_Callback(hObject, eventdata, handles)
% hObject    handle to browse_dce (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Update handles structure
guidata(hObject, handles);

[filename, pathname, filterindex] = uigetfile( ...
    {  '*.nii','Nifti Files (*.nii)'; ...
	'*2dseq','Bruker Files (2dseq)'; ...
    '*.hdr;*.img','Analyze Files (*.hdr, *.img)';...
    '*.*',  'All Files (*.*)'}, ...
    'Choose DCE-MRI file'); %#ok<NASGU>
if isequal(filename,0)
    %disp('User selected Cancel')
else
    %disp(['User selected ', fullfile(pathname, filename)])
    
    % Combine path and filename together
    fullpath = strcat(pathname,filename);

    set(handles.dce_path,'String',fullpath);
end
guidata(hObject, handles);


function dce_path_Callback(hObject, eventdata, handles)
% hObject    handle to dce_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dce_path as text
%        str2double(get(hObject,'String')) returns contents of dce_path as a double


% --- Executes during object creation, after setting all properties.
function dce_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dce_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_aif.
function browse_aif_Callback(hObject, eventdata, handles)
% hObject    handle to browse_aif (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);

[filename, pathname, filterindex] = uigetfile( ...
    {  '*.nii','Nifti Files (*.nii)'; ...
	'*2dseq','Bruker Files (2dseq)'; ...
    '*.hdr;*.img','Analyze Files (*.hdr, *.img)';...
    '*.*',  'All Files (*.*)'}, ...
    'Choose T1 map of AIF or Ref Region'); %#ok<NASGU>
if isequal(filename,0)
    %disp('User selected Cancel')
else
    %disp(['User selected ', fullfile(pathname, filename)])
    
    % Combine path and filename together
    fullpath = strcat(pathname,filename);

    set(handles.t1_aif_path,'String',fullpath);
end
guidata(hObject, handles);


function t1_aif_path_Callback(hObject, eventdata, handles)
% hObject    handle to t1_aif_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t1_aif_path as text
%        str2double(get(hObject,'String')) returns contents of t1_aif_path as a double


% --- Executes during object creation, after setting all properties.
function t1_aif_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t1_aif_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_roi.
function browse_roi_Callback(hObject, eventdata, handles)
% hObject    handle to browse_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);

[filename, pathname, filterindex] = uigetfile( ...
    {  '*.nii','Nifti Files (*.nii)'; ...
	'*2dseq','Bruker Files (2dseq)'; ...
    '*.hdr;*.img','Analyze Files (*.hdr, *.img)';...
    '*.*',  'All Files (*.*)'}, ...
    'Choose T1 map of ROI'); %#ok<NASGU>
if isequal(filename,0)
    %disp('User selected Cancel')
else
    %disp(['User selected ', fullfile(pathname, filename)])
    
    % Combine path and filename together
    fullpath = strcat(pathname,filename);

    set(handles.t1_roi_path,'String',fullpath);
end
guidata(hObject, handles);


function t1_roi_path_Callback(hObject, eventdata, handles)
% hObject    handle to t1_roi_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t1_roi_path as text
%        str2double(get(hObject,'String')) returns contents of t1_roi_path as a double


% --- Executes during object creation, after setting all properties.
function t1_roi_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t1_roi_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_noise.
function browse_noise_Callback(hObject, eventdata, handles)
% hObject    handle to browse_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);

[filename, pathname, filterindex] = uigetfile( ...
    {  '*.nii','Nifti Files (*.nii)'; ...
	'*2dseq','Bruker Files (2dseq)'; ...
    '*.hdr;*.img','Analyze Files (*.hdr, *.img)';...
    '*.*',  'All Files (*.*)'}, ...
    'Choose Noise Region'); %#ok<NASGU>
if isequal(filename,0)
    %disp('User selected Cancel')
else
    %disp(['User selected ', fullfile(pathname, filename)])
    
    % Combine path and filename together
    fullpath = strcat(pathname,filename);

    set(handles.noise_path,'String',fullpath);
end
guidata(hObject, handles);


function noise_path_Callback(hObject, eventdata, handles)
% hObject    handle to noise_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of noise_path as text
%        str2double(get(hObject,'String')) returns contents of noise_path as a double


% --- Executes during object creation, after setting all properties.
function noise_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noise_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in run_a.
function run_a_Callback(hObject, eventdata, handles)
% hObject    handle to run_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('User selected Run A')
dce_path = get(handles.dce_path,'String');
t1_aif_path = get(handles.t1_aif_path,'String');
t1_roi_path = get(handles.t1_roi_path,'String');
noise_path = get(handles.noise_path,'String');
tr = str2num(get(handles.tr, 'String')); %#ok<ST2NM>
fa = str2num(get(handles.fa, 'String')); %#ok<ST2NM>
time_resolution = str2num(get(handles.time_resolution, 'String')); %#ok<ST2NM>
time_resolution = time_resolution/60; %convert to minutes
hematocrit = str2num(get(handles.hematocrit, 'String')); %#ok<ST2NM>
snr_filter = str2num(get(handles.snr_filter, 'String')); %#ok<ST2NM>
relaxivity = str2num(get(handles.relaxivity, 'String')); %#ok<ST2NM>
injection_time = str2num(get(handles.injection_time, 'String')); %#ok<ST2NM>
water_fraction = str2num(get(handles.water_fraction, 'String')); %#ok<ST2NM>

saved_results = A_make_R1maps_func(dce_path,t1_aif_path,t1_roi_path,noise_path,tr,fa,time_resolution,hematocrit,snr_filter,relaxivity,injection_time,water_fraction);
% saved_results = 'aaa';
set(handles.results_a_path,'String',saved_results);

function tr_Callback(hObject, eventdata, handles)
% hObject    handle to tr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tr as text
%        str2double(get(hObject,'String')) returns contents of tr as a double


% --- Executes during object creation, after setting all properties.
function tr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hematocrit_Callback(hObject, eventdata, handles)
% hObject    handle to hematocrit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hematocrit as text
%        str2double(get(hObject,'String')) returns contents of hematocrit as a double


% --- Executes during object creation, after setting all properties.
function hematocrit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hematocrit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function water_fraction_Callback(hObject, eventdata, handles)
% hObject    handle to water_fraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of water_fraction as text
%        str2double(get(hObject,'String')) returns contents of water_fraction as a double


% --- Executes during object creation, after setting all properties.
function water_fraction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to water_fraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fa_Callback(hObject, eventdata, handles)
% hObject    handle to fa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fa as text
%        str2double(get(hObject,'String')) returns contents of fa as a double


% --- Executes during object creation, after setting all properties.
function fa_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function snr_filter_Callback(hObject, eventdata, handles)
% hObject    handle to snr_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of snr_filter as text
%        str2double(get(hObject,'String')) returns contents of snr_filter as a double


% --- Executes during object creation, after setting all properties.
function snr_filter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to snr_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function injection_time_Callback(hObject, eventdata, handles)
% hObject    handle to injection_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of injection_time as text
%        str2double(get(hObject,'String')) returns contents of injection_time as a double


% --- Executes during object creation, after setting all properties.
function injection_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to injection_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function time_resolution_Callback(hObject, eventdata, handles)
% hObject    handle to time_resolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_resolution as text
%        str2double(get(hObject,'String')) returns contents of time_resolution as a double


% --- Executes during object creation, after setting all properties.
function time_resolution_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_resolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function relaxivity_Callback(hObject, eventdata, handles)
% hObject    handle to relaxivity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of relaxivity as text
%        str2double(get(hObject,'String')) returns contents of relaxivity as a double


% --- Executes during object creation, after setting all properties.
function relaxivity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to relaxivity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_results_a.
function browse_results_a_Callback(hObject, eventdata, handles)
% hObject    handle to browse_results_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);

[filename, pathname, filterindex] = uigetfile( ...
    {  '*.mat','Matlab Worksapce Files (*.mat)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Choose Results from Section A'); %#ok<NASGU>
if isequal(filename,0)
    %disp('User selected Cancel')
else
    %disp(['User selected ', fullfile(pathname, filename)])
    
    % Combine path and filename together
    fullpath = strcat(pathname,filename);

    set(handles.results_a_path,'String',fullpath);
end
guidata(hObject, handles);


function results_a_path_Callback(hObject, eventdata, handles)
% hObject    handle to results_a_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of results_a_path as text
%        str2double(get(hObject,'String')) returns contents of results_a_path as a double


% --- Executes during object creation, after setting all properties.
function results_a_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to results_a_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in average_aif.
function average_aif_Callback(hObject, eventdata, handles)
% hObject    handle to average_aif (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of average_aif


% --- Executes on button press in fit_aif.
function fit_aif_Callback(hObject, eventdata, handles)
% hObject    handle to fit_aif (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fit_aif



function start_time_Callback(hObject, eventdata, handles)
% hObject    handle to start_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of start_time as text
%        str2double(get(hObject,'String')) returns contents of start_time as a double


% --- Executes during object creation, after setting all properties.
function start_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function end_time_Callback(hObject, eventdata, handles)
% hObject    handle to end_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of end_time as text
%        str2double(get(hObject,'String')) returns contents of end_time as a double


% --- Executes during object creation, after setting all properties.
function end_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to end_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function start_injection_Callback(hObject, eventdata, handles)
% hObject    handle to start_injection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of start_injection as text
%        str2double(get(hObject,'String')) returns contents of start_injection as a double


% --- Executes during object creation, after setting all properties.
function start_injection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start_injection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function end_injection_Callback(hObject, eventdata, handles)
% hObject    handle to end_injection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of end_injection as text
%        str2double(get(hObject,'String')) returns contents of end_injection as a double


% --- Executes during object creation, after setting all properties.
function end_injection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to end_injection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in run_b.
function run_b_Callback(hObject, eventdata, handles)
% hObject    handle to run_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('User selected Run B')
results_a_path = get(handles.results_a_path,'String');

start_time = str2num(get(handles.start_time, 'String')); %#ok<ST2NM>
end_time = str2num(get(handles.end_time, 'String')); %#ok<ST2NM>
start_injection = str2num(get(handles.start_injection, 'String')); %#ok<ST2NM>
end_injection = str2num(get(handles.end_injection, 'String')); %#ok<ST2NM>

fit_aif = get(handles.fit_aif, 'Value'); 
average_aif = get(handles.average_aif, 'Value'); 

saved_results = B_AIF_fitting_func(results_a_path,start_time,end_time,start_injection,end_injection,fit_aif,average_aif);
% saved_results = 'aaaa';
set(handles.results_b_path,'String',saved_results);

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


% --- Executes on button press in neuroecon.
function neuroecon_Callback(hObject, eventdata, handles)
% hObject    handle to neuroecon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of neuroecon



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



function xy_smooth_size_Callback(hObject, eventdata, handles)
% hObject    handle to xy_smooth_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xy_smooth_size as text
%        str2double(get(hObject,'String')) returns contents of xy_smooth_size as a double


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



function time_smoothing_window_Callback(hObject, eventdata, handles)
% hObject    handle to time_smoothing_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_smoothing_window as text
%        str2double(get(hObject,'String')) returns contents of time_smoothing_window as a double


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


% --- Executes on button press in browse_results_d.
function browse_results_d_Callback(hObject, eventdata, handles)
% hObject    handle to browse_results_d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);

[filename, pathname, filterindex] = uigetfile( ...
    {  '*.mat','Matlab Worksapce Files (*.mat)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Choose Results from Section D'); %#ok<NASGU>
if isequal(filename,0)
    %disp('User selected Cancel')
else
    %disp(['User selected ', fullfile(pathname, filename)])
    
    % Combine path and filename together
    fullpath = strcat(pathname,filename);

    set(handles.results_d_path,'String',fullpath);
end
guidata(hObject, handles);


function results_d_path_Callback(hObject, eventdata, handles)
% hObject    handle to results_d_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of results_d_path as text
%        str2double(get(hObject,'String')) returns contents of results_d_path as a double


% --- Executes during object creation, after setting all properties.
function results_d_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to results_d_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_background_image.
function browse_background_image_Callback(hObject, eventdata, handles)
% hObject    handle to browse_background_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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
guidata(hObject, handles);


function background_image_path_Callback(hObject, eventdata, handles)
% hObject    handle to background_image_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of background_image_path as text
%        str2double(get(hObject,'String')) returns contents of background_image_path as a double


% --- Executes during object creation, after setting all properties.
function background_image_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to background_image_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in run_d.
function run_d_Callback(hObject, eventdata, handles)
% hObject    handle to run_d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('User selected Run D')
results_b_path = get(handles.results_b_path,'String');

dce_model = get(get(handles.dce_model,'SelectedObject'),'Tag');
time_smoothing = get(get(handles.time_smoothing,'SelectedObject'),'Tag');

time_smoothing_window = str2num(get(handles.time_smoothing_window, 'String')); %#ok<ST2NM>
xy_smooth_size = str2num(get(handles.xy_smooth_size, 'String')); %#ok<ST2NM>
number_cpus = str2num(get(handles.number_cpus, 'String')); %#ok<ST2NM>

neuroecon = get(handles.neuroecon, 'Value'); 


saved_results = D_fit_voxels_func(results_b_path,dce_model,time_smoothing,time_smoothing_window,xy_smooth_size,number_cpus,neuroecon);
% saved_results = 'aaaa';
set(handles.results_d_path,'String',saved_results);

% --- Executes on button press in run_e.
function run_e_Callback(hObject, eventdata, handles)
% hObject    handle to run_e (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

results_d_path = get(handles.results_d_path,'String');
background_image_path = get(handles.background_image_path,'String');


compare_fits(results_d_path,background_image_path);

