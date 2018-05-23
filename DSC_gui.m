function varargout = DSC_gui(varargin)
% DSC_GUI MATLAB code for DSC_gui.fig
%      DSC_GUI, by itself, creates a new DSC_GUI or raises the existing
%      singleton*.
%
%      H = DSC_GUI returns the handle to a new DSC_GUI or the handle to
%      the existing singleton*.
%
%      DSC_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DSC_GUI.M with the given input arguments.
%
%      DSC_GUI('Property','Value',...) creates a new DSC_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DSC_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DSC_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DSC_gui

% Last Modified by GUIDE v2.5 17-May-2018 14:26:40

% REVISION HISTORY: 
% 04/26/2015:AIF and NOISE handling optins have been added.  
% 04/27/2015: Fitting will be dialed in for both humans and  mouse.  
% 04/27/2015: Set defaults: Auto AIF; Auto noise; species, 'mouse'. 
% Tested GUI for both human and mouse, slice images, using auto AIF/Auto
% noise. Still have to debug manual AIF/ manual noise for both species
% using slices. Also, stacks remain to be check for all conditions.  

% Dependent functions: 
% DSC_signal2concentration: 
% AIF_auto_cluster:
% AIF_manual: 
% gfun: 
% DSC_convolution_sSVD
% parse_preferene_file (part of Rocketship package). 

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DSC_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @DSC_gui_OutputFcn, ...
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


% --- Executes just before DSC_gui is made visible.
function DSC_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DSC_gui (see VARARGIN)

% Choose default command line output for DSC_gui
handles.output = hObject;

%ADDITIONAL BOXES IN THE OUTPUT STRUCTURE: 
handles.dsc_image = []; 
handles.AIF_roi = [];
handles.import_aif = [];
handles.import_store = [];
handles.noise_roi = []; 
handles.TE = [];
handles.TR = []; 
handles.rho = []; 
handles.Psvd = []; 
handles.species = 'mouse'; 
handles.r2_star = []; 
handles.AIF_type = 0; 
handles.noise_type = 0;
handles.Fitting_function = 0;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DSC_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DSC_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in add_images.
function add_images_Callback(hObject, eventdata, handles)
% hObject    handle to add_images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% When this button gets pressed, open a dialog window to select a file: 
[filename, pathname, filterspec] = uigetfile('*.nii','Nifti Files (*.nii)'); 

if isequal(filename,0)
    %disp('User selected Cancel')
    fullpath = '';
else
    %disp(['User selected ', fullfile(pathname, filename)])
    
    % Combine path and filename together
    fullpath = strcat(pathname,filename);
    path_out = fullpath; 
    
    if ischar(fullpath)
        fullpath = {fullpath};
    end
    % if isempty(list)
    %     list = {''};
    % end
    
    %filename = filename';
    fullpath = fullpath';
 
    
    if numel(fullpath) > 1
        visualpath = [fullpath{1} ' -> ' num2str(numel(fullpath)) ' files'];
    else
        visualpath = fullpath{1};
    end
    
    set(handles.dsc_image_path,'String',visualpath);
    handles.dsc_image = path_out; 
  
end

guidata(hObject, handles);


% --- Executes on button press in remoive_images.
function remoive_images_Callback(hObject, eventdata, handles)
% hObject    handle to remoive_images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.dsc_image_path = []; 
handles.dsc_image = []; 
guidata(hObject, handles);


function TE_input_Callback(hObject, eventdata, handles)
% hObject    handle to TE_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TE_input as text
%        str2double(get(hObject,'String')) returns contents of TE_input as a double
handles.TE = str2double(get(hObject,'String')); 
uiremember; 
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function TE_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TE_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
uirestore; 


function TR_input_Callback(hObject, eventdata, handles)
% hObject    handle to TR_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TR_input as text
%        str2double(get(hObject,'String')) returns contents of TR_input as a double
handles.TR = str2double(get(hObject,'String')); 
uiremember; 
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function TR_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TR_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
uirestore; 


function rho_input_Callback(hObject, eventdata, handles)
% hObject    handle to rho_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rho_input as text
%        str2double(get(hObject,'String')) returns contents of rho_input as a double
handles.rho = str2double(get(hObject,'String')); 
uiremember; 
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function rho_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rho_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
uirestore; 

function r2_star_input_Callback(hObject, eventdata, handles)
% hObject    handle to r2_star_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of r2_star_input as text
%        str2double(get(hObject,'String')) returns contents of r2_star_input as a double
handles.r2_star = str2double(get(hObject,'String')); 
uiremember; 
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function r2_star_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r2_star_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
uirestore;



function Psvd_input_Callback(hObject, eventdata, handles)
% hObject    handle to Psvd_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Psvd_input as text
%        str2double(get(hObject,'String')) returns contents of Psvd_input as a double
handles.Psvd = str2double(get(hObject,'String')); 
uiremember; 
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Psvd_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Psvd_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
uirestore;


function AIF_path_Callback(hObject, eventdata, handles)
% hObject    handle to AIF_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AIF_path as text
%        str2double(get(hObject,'String')) returns contents of AIF_path as a double


% --- Executes during object creation, after setting all properties.
function AIF_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AIF_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
% --- Executes on button press in add_AIF_roi.
function add_AIF_roi_Callback(hObject, eventdata, handles)
% hObject    handle to add_AIF_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% When this button gets pressed, open a dialog window to select a file: 
[filename, pathname, filterspec] = uigetfile('*.nii','Nifti Files (*.nii)'); 

if isequal(filename,0)
    %disp('User selected Cancel')
    fullpath = '';
else
    %disp(['User selected ', fullfile(pathname, filename)])
    
    % Combine path and filename together
    fullpath = strcat(pathname,filename);
    path_out = fullpath; 
    
    if ischar(fullpath)
        fullpath = {fullpath};
    end
    % if isempty(list)
    %     list = {''};
    % end
    
    %filename = filename';
    fullpath = fullpath';
    
    if numel(fullpath) > 1
        visualpath = [fullpath{1} ' -> ' num2str(numel(fullpath)) ' files'];
    else
        visualpath = fullpath{1};
    end
    
    set(handles.AIF_path,'String',visualpath);
    handles.AIF_roi = path_out; 
   
end
guidata(hObject, handles);


function import_path_Callback(hObject, eventdata, handles)
% hObject    handle to import_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of import_path as text
%        str2double(get(hObject,'String')) returns contents of import_path as a double


% --- Executes during object creation, after setting all properties.
function import_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to import_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in add_import_aif.
function add_import_aif_Callback(hObject, eventdata, handles)
% hObject    handle to add_import_aif (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname, filterspec] = uigetfile('*.nii','Nifti Files (*.nii)'); 

if isequal(filename,0)
    %disp('User selected Cancel')
    fullpath = '';
else
    %disp(['User selected ', fullfile(pathname, filename)])
    
    % Combine path and filename together
    fullpath = strcat(pathname,filename);
    path_out = fullpath; 
    
    if ischar(fullpath)
        fullpath = {fullpath};
    end
    % if isempty(list)
    %     list = {''};
    % end
    
    %filename = filename';
    fullpath = fullpath';
    
    if numel(fullpath) > 1
        visualpath = [fullpath{1} ' -> ' num2str(numel(fullpath)) ' files'];
    else
        visualpath = fullpath{1};
    end
    
    set(handles.import_path,'String',visualpath);
    handles.AIF_roi = path_out; 
   
end
guidata(hObject, handles);

% --- Executes on button press in Import_stored.
function Import_stored_Callback(hObject, eventdata, handles)
% hObject    handle to Import_stored (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Import_stored
[filename, pathname, filterspec] = uigetfile('*.nii','Nifti Files (*.nii)'); 
guidata(hObject, handles);

function dsc_image_path_Callback(hObject, eventdata, handles)
% hObject    handle to dsc_image_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dsc_image_path as text
%        str2double(get(hObject,'String')) returns contents of dsc_image_path as a double


% --- Executes during object creation, after setting all properties.
function dsc_image_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dsc_image_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when selected object is changed in species_selection.
function species_selection_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in species_selection 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
% Settingn the initial case: 
set(handles.human_button,'Value', 0); 
set(handles.mouse_button,'Value',1); 
switch get(eventdata.NewValue, 'Tag') 
    case 'human_button' 
        handles.species = 'human'; 
        set(handles.human_button,'Value', 1); 
        set(handles.mouse_button,'Value',0);
    case 'mouse_button' 
        handles.species = 'mouse'; 
        set(handles.human_button,'Value', 0); 
        set(handles.mouse_button,'Value',1); 
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function species_selection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to species_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function mouse_button_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mouse_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function human_button_CreateFcn(hObject, eventdata, handles)
% hObject    handle to human_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


%SETUP IS DONE. NOW WE SEND THE COLLECTED INORMATION TO SUBSEQUENT
%FUNCTIONS

% --- Executes on button press in Start_button.
function Start_button_Callback(hObject, eventdata, handles)
% hObject    handle to Start_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Unboxing the variables from the handles structure: 
TE = str2double(get(handles.TE_input,'String'))/1000.0; %conversion from miliseconds to seconds 12/15/17
TR = str2double(get(handles.TR_input,'String'));%/1000.0;
deltaT = TR / 60;  % converstion to minutes. 
r2_star = str2double(get(handles.r2_star_input,'String')); 
species = handles.species; 
Psvd = str2double(get(handles.Psvd_input,'String')); 
rho =  str2double(get(handles.rho_input,'String'));
AIF_type = handles.AIF_type; 
Fitting_function = handles.Fitting_function;
noise_type = handles.noise_type; 
details = whos('r2_star'); 

%Next WE ARE GOING TO LOAD THE IMAGE FILE: 
disp(handles.dsc_image); 
disp(handles.species); 
disp(handles.r2_star); 
image_array = load_nii(handles.dsc_image); 
image_array = image_array.img; 
%Next we load the nosie roi array: 
if noise_type == 0 
    roi_array = []; 
elseif noise_type == 1 
    roi_array = load_nii(handles.noise_roi); 
    roi_array = roi_array.img; 
end


%getting the path from the image file: 
[ ~ , image_path] = fileparts(handles.dsc_image); 

[concentration_array, base_concentration_array, time_vect, base_time_vect,whole_time_vect, bolus_time] = DSC_signal2concentration(image_array,TE,TR,r2_star,species,image_path,noise_type, ...
   roi_array ); 
if AIF_type == 0 %AIF Auto
    [meanAIF, meanSignal] = AIF_auto_cluster(concentration_array, image_array, time_vect, TR,species);
    baseline = 0; %develop a way to calculate the baseline from the auto clustered AIF
    baseline_array = zeros(numel(base_time_vect));
elseif AIF_type ==1 %AIF User Selected
    AIF_mask = load_nii(handles.AIF_roi); 
    AIF_mask = AIF_mask.img; 
    [meanAIF, meanSignal, baseline,baseline_array] = AIF_manual_noreshape(image_array,concentration_array,AIF_mask, base_time_vect,base_concentration_array);

elseif AIF_type==2 %AIF Import
    load(handles.import_AIF)
    [meanAIF_adjusted, time_vect, concentration_array] = import_AIF(meanAIF, bolus_time, time_vect, concentration_array, r2_star, TE);
    meanAIF = meanAIF_adjusted;
    
    display('time_vect_num')
    numel(time_vect)
    display('meanAIF_num')
    numel(meanAIF)
    
elseif AIF_type==3 %AIF Use Previous
    load('previous_data.mat');
    [meanAIF_adjusted, time_vect, concentration_array] = previous_AIF(meanAIF,meanSignal,bolus_time, time_vect,concentration_array);
    meanAIF = meanAIF_adjusted;
    
    display('time_vect_num')
    numel(time_vect)
    display('meanAIF_num')
    numel(meanAIF)
end

%create a mean AIF spaning the whole scan time 
AIF_whole = cat(1,baseline_array, meanAIF);

%now run the selected fitting function
if Fitting_function == 0 %forced linear biexponential (uses local max) %the upslope is fitted to
    
    Cp = cat(1,baseline_array,meanAIF);
    step = [(bolus_time) (bolus_time + numel(time_vect))];
    T1 = whole_time_vect;
    xdata = struct('Cp',Cp,'baseline', baseline, 'timer', T1, 'step', step, 'bolus_time', bolus_time);
    verbose = -1; %this prevents any internal verbose function from running change to 1 to run verbose

    [Cp, x, xdata, rsqurare] = Single_Forced_linear_AIFbiexpfithelplocal(xdata,verbose);
    Ct = Cp;
    Ct(1:bolus_time - 1) = [];
    Ct = Ct';
    
elseif Fitting_function == 1 %biexponential (uses absolute max)
    Cp = cat(1,baseline_array,meanAIF);
    step = [bolus_time (bolus_time + numel(time_vect))];
    T1 = whole_time_vect;
    xdata = struct('Cp',Cp,'baseline', baseline, 'timer', T1, 'step', step, 'bolus_time', bolus_time);
    verbose = -1; %this prevents any internal verbose function from running change to 1 to run verbose

    [Cp, x, xdata, rsqurare] = AIFbiexpfithelp(xdata,verbose);
    Ct = Cp;
    Ct(1:bolus_time - 1) = [];
    Ct = Ct';
    
elseif Fitting_function == 2 %biexponential (uses local max) 
    Cp = cat(1,baseline_array,meanAIF);
    step = [bolus_time (bolus_time + numel(time_vect))];
    T1 = whole_time_vect;
    xdata = struct('Cp',Cp,'baseline', baseline, 'timer', T1, 'step', step, 'bolus_time', bolus_time);
    verbose = -1; %this prevents any internal verbose function from running change to 1 to run verbose

    [Cp, x, xdata, rsqurare] = AIFbiexpfithelplocal(xdata,verbose);
    Ct = Cp;
    Ct(1:bolus_time - 1) = [];
    Ct = Ct';
    
elseif Fitting_function == 3 %gamma-variant
    % Now we fit the AIF with a SCR model: 

    %assigning the gamma variate function, gfun, to be our desired fitting
    %function: 
    Ct = fitting_gamma_variant(meanAIF,species, time_vect);
    Cp = cat(1,baseline_array, Ct); %Cp is created for plotting purposes only. Ct is analyzed for CBF, CBV...
    
elseif Fitting_function == 4 %raw data 
    Ct = meanAIF;
    Cp = cat(1,baseline_array, Ct);

%{
elseif Fitting_function == 42 %copy_of_upslopewith peak based decision making (not currently an option)
    Cp = cat(1,baseline_array,meanAIF);
    step = [bolus_time (bolus_time + numel(time_vect))];
    T1 = whole_time_vect;
    xdata = struct('Cp',Cp,'baseline', baseline, 'timer', T1, 'step', step, 'bolus_time', bolus_time);
    verbose = -1; %this prevents any internal verbose function from running change to 1 to run verbose

    [Cp, x, xdata, rsqurare] = Forced_linear_AIFbiexpfithelplocal(xdata,verbose);
    Ct = Cp;
    Ct(1:bolus_time - 1) = [];
    [~,max_indexCt] = max(Ct); 
    
    if numel(findpeaks(AIF_whole(bolus_time:((bolus_time + 1) + max_indexCt)))) == 1 %check for additional local maxima in upslope, if none
        %are present use the exact upslope vales for the fitted funciton
        Ct(1:max_indexCt) = meanAIF(1: max_indexCt);
        Ct = Ct';
        base = baseline * ones(numel(baseline_array),1);
        Cp = [base; Ct];
    else %if additional local maxima or 'bumps' in the upslope are present, use the 2 point forced linear fit
        Ct = Ct';
    end
%}
elseif Fitting_function == 5 %copy_of_upslope with peak based decision making
    Cp = cat(1,baseline_array,meanAIF);
    step = [bolus_time (bolus_time + numel(time_vect))];
    T1 = whole_time_vect;
    xdata = struct('Cp',Cp,'baseline', baseline, 'timer', T1, 'step', step, 'bolus_time', bolus_time);
    verbose = -1; %this prevents any internal verbose function from running change to 1 to run verbose

    [Cp, x, xdata, rsqurare] = AIFbiexpfithelplocal(xdata,verbose);
    Ct = Cp;
    Ct(1:bolus_time -1) = [];
    
   %calculate the most likely first peak by finding the first local maxima
    %as was done in the local fitting help function
    %peak is the first within 3% of the max value of the whole data set
    [local_maxima, maxima_indexes] = findpeaks(Ct);
    maxima_iterator = 1;
    while(local_maxima(maxima_iterator) < (0.97 * max(local_maxima)))
        maxima_iterator = maxima_iterator + 1;
    end

    max_indexCt = maxima_indexes(maxima_iterator);
    
    if numel(findpeaks(AIF_whole(bolus_time:(bolus_time + max_indexCt)))) == 1 %check for additional local maxima in upslope, if none
        %are present use the exact upslope vales for the fitted funciton
        %[bolus time, max_index + 1] is the range that is examined for the
        %additional maxima
        Ct(1:max_indexCt) = meanAIF(1:max_indexCt);
        Ct = Ct';
        base = baseline * ones(numel(baseline_array),1);
        Cp = [base; Ct];
    else %if additional local maxima or 'bumps' in the upslope are present, use a smooth fit
        Ct = Ct';
    end
    
elseif Fitting_function == 6 %upslope copy biexponetial
    
    Cp = cat(1,baseline_array,meanAIF);
    step = [bolus_time (bolus_time + numel(time_vect))];
    T1 = whole_time_vect;
    xdata = struct('Cp',Cp,'baseline', baseline, 'timer', T1, 'step', step, 'bolus_time', bolus_time);
    verbose = -1; %this prevents any internal verbose function from running change to 1 to run verbose

    [Cp, x, xdata, rsqurare] = AIFbiexpfithelplocal(xdata,verbose);
    Ct = Cp;
    Ct(1:bolus_time - 1) = [];
    
    %calculate the most likely first peak by finding the first local maxima
    %as was done in the local fitting help function
    %peak is the first within 3% of the max value of the whole data set
    [local_maxima, maxima_indexes] = findpeaks(Ct);
    maxima_iterator = 1;
    while(local_maxima(maxima_iterator) < (0.97 * max(local_maxima)))
        maxima_iterator = maxima_iterator + 1;
    end

    max_indexCt = maxima_indexes(maxima_iterator);

    max_indexCt = maxima_indexes(maxima_iterator);
    
    %exact upslope vales for the fitted funciton
    Ct(1:max_indexCt) = meanAIF(1:max_indexCt);
    Ct = Ct';
    base = baseline * ones(numel(baseline_array),1);
    Cp = [base; Ct];
    
elseif Fitting_function == 7 %forced bilinear biexponential decay
    
    %two linear function
    Cp = cat(1,baseline_array,meanAIF);
    step = [(bolus_time) (bolus_time + numel(time_vect))];
    T1 = whole_time_vect;
    xdata = struct('Cp',Cp,'baseline', baseline, 'timer', T1, 'step', step, 'bolus_time', bolus_time);
    verbose = -1; %this prevents any internal verbose function from running change to 1 to run verbose

    [Cp, x, xdata, rsqurare] = Forced_linear_AIFbiexpfithelplocal(xdata,verbose);
    Ct = Cp;
    Ct(1:bolus_time - 1) = [];
    Ct = Ct';
    
end

%save this runs meanAIF, bolus time, and importedAIF
save('previous_data.mat', 'meanAIF','meanSignal','bolus_time','baseline');

%plot the meanAIF and fitted function over the entire scan
figure;
plot(whole_time_vect,Cp,'-o', whole_time_vect, AIF_whole);
title('AIF and Fitted AIF Over time');
xlabel('Time (min)');                                                                             
ylabel('Concentration (mM)');
legend('Fitted AIF', 'Raw AIF');

%plot the mean AIF and fitted function over the time period to be analyzed
figure; 
time_vect_sec = time_vect * 60; %get time in seconds to plot
time_vect_sec = time_vect_sec + (bolus_time);
%base_time_vect_sec = base_time_vect * 60;
plot(time_vect_sec,Ct,'-o',time_vect_sec,meanAIF); 
title('AIF and Fitted AIF Over time');
xlabel('Time (s)');                                                                             
ylabel('Concentration (mM)');
legend('Fitted AIF', 'Raw AIF');

%plot the mean AIF signal
figure; 
time_vect_sec2 = 0 : length(meanSignal) -1; 
time_vect_sec2 = time_vect_sec2 * TR;
plot(time_vect_sec2, meanSignal, 'b');
title('AIF Mean Signal from Start to End of Scan')
xlabel('Time (s)')
ylabel('Signal Intensity (au)')

Kh = 0.71; 
% method = 1; 
[CBF, CBV, MTT] = DSC_convolution_sSVD(concentration_array,Ct,deltaT,Kh,rho,Psvd,1,image_path); 
disp('new')
% guidata(hObject, handles);



% --- Executes on button press in Cancel_button.
function Cancel_button_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in add_noise_roi.
function add_noise_roi_Callback(hObject, eventdata, handles)
% hObject    handle to add_noise_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname, filterspec] = uigetfile('*.nii','Nifti Files (*.nii)'); 

if isequal(filename,0)
    %disp('User selected Cancel')
    fullpath = '';
else
    %disp(['User selected ', fullfile(pathname, filename)])
    
    % Combine path and filename together
    fullpath = strcat(pathname,filename);
    path_out = fullpath; 
    
    if ischar(fullpath)
        fullpath = {fullpath};
    end
    % if isempty(list)
    %     list = {''};
    % end
    
    %filename = filename';
    fullpath = fullpath';
    
    if numel(fullpath) > 1
        visualpath = [fullpath{1} ' -> ' num2str(numel(fullpath)) ' files'];
    else
        visualpath = fullpath{1};
    end
    
    set(handles.noise_file,'String',visualpath);
    handles.noise_roi = path_out; 
   
end
guidata(hObject, handles);



function noise_file_Callback(hObject, eventdata, handles)
% hObject    handle to noise_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of noise_file as text
%        str2double(get(hObject,'String')) returns contents of noise_file as a double
% Do not enable unless button press...
set(handles.noise_file,'enable', 'off'); 
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function noise_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noise_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function user_Noise_CreateFcn(hObject, eventdata, handles)
% hObject    handle to user_Noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function auto_Noise_CreateFcn(hObject, eventdata, handles)
% hObject    handle to user_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when selected object is changed in noise_selection.
function noise_selection_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in noise_selection 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
set(handles.user_Noise,'Value',0); 
set(handles.auto_Noise,'Value',1); 
switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'auto_Noise'
        handles.noise_type =  0; 
        set(handles.user_Noise,'Value',0); 
        set(handles.auto_Noise,'Value',1);
        set(handles.add_noise_roi,'enable', 'off'); 
        set(handles.noise_file,'enable', 'off'); 
    case 'user_Noise'
        handles.noise_type = 1; 
        set(handles.user_Noise,'Value',1); 
        set(handles.auto_Noise,'Value',0);
        set(handles.add_noise_roi,'enable', 'on'); 
        set(handles.noise_file,'enable', 'on');  
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function auto_cluster_CreateFcn(hObject, eventdata, handles)
% hObject    handle to auto_cluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function user_aif_CreateFcn(hObject, eventdata, handles)
% hObject    handle to user_aif (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function import_aif1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to import_aif1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function Import_stored_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Import_stored (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when selected object is changed in aif_selection.
function aif_selection_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in aif_selection 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
set(handles.auto_cluster,'Value',1); 
set(handles.user_aif,'Value',0); 
set(handles.import_aif1,'Value',0); 
set(handles.add_AIF_roi,'enable', 'off'); 
set(handles.add_import_aif,'enable', 'off');
set(handles.AIF_path,'enable', 'off'); 
set(handles.import_path,'enable', 'off');
set(handles.import_store,'Value',0);

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'auto_cluster'
        handles.AIF_type =  0; 
        set(handles.add_AIF_roi,'enable', 'off'); 
        set(handles.AIF_path,'enable', 'off'); 
        set(handles.import_store,'enable', 'off');
        set(handles.auto_cluster,'Value',1); 
        set(handles.user_aif,'Value',0);
        set(handles.import_aif1,'Value',0); 
        set(handles.add_import_aif,'enable', 'off');
        set(handles.import_path,'enable', 'off');
        set(handles.import_store,'Value',0);
    case 'user_aif'
       handles.AIF_type = 1; 
       set(handles.add_AIF_roi,'enable', 'on'); 
       set(handles.AIF_path,'enable', 'on');
       set(handles.auto_cluster,'Value',0); 
       set(handles.user_aif,'Value',1);
       set(handles.import_aif1,'Value',0);
       set(handles.add_import_aif,'enable', 'off');
       set(handles.import_path,'enable', 'off');
       set(handles.import_store,'Value',0);
       set(handles.import_store,'enable', 'off');
    case 'import_aif1'
       handles.AIF_type = 2; 
       set(handles.add_AIF_roi,'enable', 'off'); 
       set(handles.AIF_path,'enable', 'off');
       set(handles.import_store,'enable', 'off');
       set(handles.import_aif1,'Value',1);
       set(handles.import_path,'enable', 'on');
       set(handles.add_import_aif,'enable', 'on');
       set(handles.auto_cluster,'Value',0); 
       set(handles.user_aif,'Value',0);
        set(handles.import_store,'Value',0);
    case 'import_stored'
       handles.AIF_type = 3; 
       set(handles.add_AIF_roi,'enable', 'off'); 
       set(handles.AIF_path,'enable', 'off');
       set(handles.import_store,'enable', 'on');
       set(handles.import_aif1,'Value',0);
       set(handles.import_store,'Value',1);
       set(handles.import_path,'enable', 'off');
       set(handles.add_import_aif,'enable', 'off');
       set(handles.auto_cluster,'Value',0); 
       set(handles.user_aif,'Value',0); 
end
guidata(hObject, handles);


% --- Executes on button press in Restart.
function Restart_Callback(hObject, eventdata, handles)
% hObject    handle to Restart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function forcelb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to forcelb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function biexp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to biexp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function biexl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to biexl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function Gammav_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gammav (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function rawdata_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rawdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes during object creation, after setting all properties.


function upslopebiexp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to upslopebiexp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function upslopebiexpliadj_CreateFcn(hObject, eventdata, handles)
% hObject    handle to upslopebiexpliadj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function forcebiexp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to forcebiexp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes when selected object is changed in Fitting_function.
function Fitting_function_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in Fitting_function 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.forcelb,'Value',1); 
set(handles.biexp,'Value',0); 
set(handles.biexl,'Value',0); 
set(handles.Gammav,'Value', 0); 
set(handles.rawdata,'Value', 0);
set(handles.upslopebiexp,'Value', 0);
set(handles.upslopebiexpliadj,'Value', 0);
set(handles.forcebiexp,'Value',0);

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'forcelb'
        handles.Fitting_function =  0; 
        set(handles.forcelb,'Value',1); 
        set(handles.biexp,'Value',0); 
        set(handles.biexl,'Value',0); 
        set(handles.Gammav,'Value', 0); 
        set(handles.rawdata,'Value', 0);
        set(handles.upslopebiexp,'Value', 0);
        set(handles.upslopebiexpliadj,'Value', 0);
        set(handles.forcebiexp,'Value',0);
    case 'biexp'
        handles.Fitting_function = 1; 
        set(handles.forcelb,'Value',0); 
        set(handles.biexp,'Value',1); 
        set(handles.biexl,'Value',0); 
        set(handles.Gammav,'Value', 0); 
        set(handles.rawdata,'Value', 0);
        set(handles.upslopebiexp,'Value', 0);
        set(handles.upslopebiexpliadj,'Value', 0);
        set(handles.forcebiexp,'Value',0);
    case 'biexl'
        handles.Fitting_function = 2; 
        set(handles.forcelb,'Value',0); 
        set(handles.biexp,'Value',0); 
        set(handles.biexl,'Value',1); 
        set(handles.Gammav,'Value', 0); 
        set(handles.rawdata,'Value', 0);
        set(handles.upslopebiexp,'Value', 0);
        set(handles.upslopebiexpliadj,'Value', 0);
        set(handles.forcebiexp,'Value',0);
    case 'Gammav'
        handles.Fitting_function = 3; 
        set(handles.forcelb,'Value',0); 
        set(handles.biexp,'Value',0); 
        set(handles.biexl,'Value',0); 
        set(handles.Gammav,'Value', 1); 
        set(handles.rawdata,'Value', 0);
        set(handles.upslopebiexp,'Value', 0);
        set(handles.upslopebiexpliadj,'Value', 0);
        set(handles.forcebiexp,'Value',0);
    case 'rawdata'
        handles.Fitting_function = 4; 
        set(handles.forcelb,'Value',0); 
        set(handles.biexp,'Value',0); 
        set(handles.biexl,'Value',0); 
        set(handles.Gammav,'Value', 0); 
        set(handles.rawdata,'Value', 1);
        set(handles.upslopebiexp,'Value', 0);
        set(handles.upslopebiexpliadj,'Value', 0);
        set(handles.forcebiexp,'Value',0);
    case 'upslopebiexp'
        handles.Fitting_function = 5; 
        set(handles.forcelb,'Value',0); 
        set(handles.biexp,'Value',0); 
        set(handles.biexl,'Value',0); 
        set(handles.Gammav,'Value', 0); 
        set(handles.rawdata,'Value', 0);
        set(handles.upslopebiexp,'Value', 1);
        set(handles.upslopebiexpliadj,'Value', 0);
        set(handles.forcebiexp,'Value',0);
    case 'upslopebiexpliadj'
        handles.Fitting_function = 6; 
        set(handles.forcelb,'Value',0); 
        set(handles.biexp,'Value',0); 
        set(handles.biexl,'Value',0); 
        set(handles.Gammav,'Value', 0); 
        set(handles.rawdata,'Value', 0);
        set(handles.upslopebiexp,'Value', 0);
        set(handles.upslopebiexpliadj,'Value', 1);
        set(handles.forcebiexp,'Value',0);
    case 'forcebiexp'
        handles.Fitting_function = 7; 
        set(handles.forcelb,'Value',0); 
        set(handles.biexp,'Value',0); 
        set(handles.biexl,'Value',0); 
        set(handles.Gammav,'Value', 0); 
        set(handles.rawdata,'Value', 0);
        set(handles.upslopebiexp,'Value', 0);
        set(handles.upslopebiexpliadj,'Value', 0);
        set(handles.forcebiexp,'Value',1);
end
guidata(hObject, handles);
