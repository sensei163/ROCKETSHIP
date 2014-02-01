function varargout = RUNB(varargin)
% RUNB MATLAB code for RUNB.fig
%      RUNB, by itself, creates a new RUNB or raises the existing
%      singleton*.
%
%      H = RUNB returns the handle to a new RUNB or the handle to
%      the existing singleton*.
%
%      RUNB('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RUNB.M with the given input arguments.
%
%      RUNB('Property','Value',...) creates a new RUNB or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RUNB_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RUNB_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RUNB

% Last Modified by GUIDE v2.5 31-Jan-2014 18:47:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RUNB_OpeningFcn, ...
                   'gui_OutputFcn',  @RUNB_OutputFcn, ...
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


% --- Executes just before RUNB is made visible.
function RUNB_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RUNB (see VARARGIN)

% Choose default command line output for RUNB
handles.output = hObject;

% Setup the name of the A results file
ARESULTS = varargin{1};
set(handles.results_a_path, 'String', ARESULTS);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RUNB wait for user response (see UIRESUME)
 uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RUNB_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} =handles.saved_results;
delete(handles.figure1);



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



function import_aif_path_Callback(hObject, eventdata, handles)
% hObject    handle to import_aif_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of import_aif_path as text
%        str2double(get(hObject,'String')) returns contents of import_aif_path as a double


% --- Executes during object creation, after setting all properties.
function import_aif_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to import_aif_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_import_aif.
function browse_import_aif_Callback(hObject, eventdata, handles)
% hObject    handle to browse_import_aif (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);

[filename, pathname, filterindex] = uigetfile( ...
    {  '*.mat','Matlab Worksapce Files (*.mat)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Choose Results with an AIF'); %#ok<NASGU>
if isequal(filename,0)
    %disp('User selected Cancel')
else
    %disp(['User selected ', fullfile(pathname, filename)])
    
    % Combine path and filename together
    fullpath = strcat(pathname,filename);

    set(handles.import_aif_path,'String',fullpath);
end
guidata(hObject, handles);


% --- Executes on selection change in aif_type.
function aif_type_Callback(hObject, eventdata, handles)
% hObject    handle to aif_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns aif_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from aif_type
uiremember;
update_importaif(handles);

function update_importaif(handles)
if get(handles.aif_type,'Value')==3
    set(handles.import_aif_path,'Enable','on');
    set(handles.browse_import_aif,'Enable','on');
else
    set(handles.import_aif_path,'Enable','off');
    set(handles.browse_import_aif,'Enable','off');
end


% --- Executes during object creation, after setting all properties.
function aif_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to aif_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
uirestore;

% --- Executes on button press in run_average_aif.
function run_average_aif_Callback(hObject, eventdata, handles)
% hObject    handle to run_average_aif (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
results_a_path = get(handles.results_a_path,'String');
[pathstr,~,~] = fileparts(results_a_path);
if ~isempty(pathstr)
    average_aifs('save_path',pathstr);
else
    average_aifs();
end



function start_time_Callback(hObject, eventdata, handles)
% hObject    handle to start_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of start_time as text
%        str2double(get(hObject,'String')) returns contents of start_time as a double
uiremember;

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
uirestore;


function end_time_Callback(hObject, eventdata, handles)
% hObject    handle to end_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of end_time as text
%        str2double(get(hObject,'String')) returns contents of end_time as a double
uiremember;

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
uirestore;


function start_injection_Callback(hObject, eventdata, handles)
% hObject    handle to start_injection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of start_injection as text
%        str2double(get(hObject,'String')) returns contents of start_injection as a double
uiremember;

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
uirestore;


function end_injection_Callback(hObject, eventdata, handles)
% hObject    handle to end_injection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of end_injection as text
%        str2double(get(hObject,'String')) returns contents of end_injection as a double
uiremember;

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
uirestore;


function time_resolution_Callback(hObject, eventdata, handles)
% hObject    handle to time_resolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_resolution as text
%        str2double(get(hObject,'String')) returns contents of time_resolution as a double
uiremember;

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
uirestore;

% --- Executes on button press in done.
function done_Callback(hObject, eventdata, handles)
% hObject    handle to done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%disp('User selected Run B');
results_a_path = get(handles.results_a_path,'String');

start_time = str2num(get(handles.start_time, 'String')); %#ok<ST2NM>
end_time = str2num(get(handles.end_time, 'String')); %#ok<ST2NM>
start_injection = str2num(get(handles.start_injection, 'String')); %#ok<ST2NM>
end_injection = str2num(get(handles.end_injection, 'String')); %#ok<ST2NM>
% fit_aif = get(handles.fit_aif, 'Value'); 
% average_aif = get(handles.average_aif, 'Value'); 
fit_aif = (1==get(handles.aif_type,'Value'));
time_resolution = str2num(get(handles.time_resolution, 'String')); %#ok<ST2NM>
time_resolution = time_resolution/60; %convert to minutes
if get(handles.aif_type,'Value')==3
    import_aif_path = get(handles.import_aif_path, 'String');
else
    import_aif_path = '';
end

saved_results = B_AIF_fitting_func(results_a_path,start_time,end_time,start_injection,end_injection,fit_aif,import_aif_path,time_resolution);
%set(handles.results_b_path,'String',saved_results);
handles.saved_results = saved_results;
guidata(hObject, handles);
uiresume(handles.figure1);

% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
% hObject    handle to cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
varargout{1} = ''; %handles.output;
delete(handles.figure1);
