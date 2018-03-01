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

% Last Modified by GUIDE v2.5 15-Feb-2018 12:17:31

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
handles.noise_roi = []; 
handles.TE = [];
handles.TR = []; 
handles.rho = []; 
handles.Psvd = []; 
handles.species = 'mouse'; 
handles.r2_star = []; 
handles.AIF_type = 0; 
handles.noise_type = 0; 

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

[concentration_array, base_concentration_array, time_vect, base_time_vect, bolus_time] = DSC_signal2concentration(image_array,TE,TR,r2_star,species,image_path,noise_type, ...
   roi_array ); 
if AIF_type == 0 
    [meanAIF, meanSignal] = AIF_auto_cluster(concentration_array, image_array, time_vect, TR,species); 
elseif AIF_type ==1 
    AIF_mask = load_nii(handles.AIF_roi); 
    AIF_mask = AIF_mask.img; 
    [meanAIF, meanSignal] = AIF_manual(image_array,concentration_array,AIF_mask);
elseif AIF_type ==2 %AIF_import file saved as a handles.import_AIF
    load('previous_meanAIF_signalIntensity_bolusTime.mat');
    [meanAIF_adjusted, time_vect, concentration_array] = previous_AIF(meanAIF,meanSignal,bolus_time, time_vect,concentration_array);
    meanAIF = meanAIF_adjusted;
    
    display('time_vect_num')
    numel(time_vect)
    display('meanAIF_num')
    numel(meanAIF)
    
    %adjust lenth of concentration array to reflect the period for which
    %we have an AIF 
    
    %{
    AIF_mask = load_nii(handles.AIF_roi);
    AIF_mask = AIF_mask.img;
    [meanAIF, meanSignal] = AIF_import(AIF_mask);
    %}
end 

%save this runs meanAIF, bolus time, and importedAIF

save('previous_meanAIF_signalIntensity_bolusTime.mat', 'meanAIF','meanSignal','bolus_time');

% Now we fit the AIF with a SCR model: 

%assigning the gamma variate function, gfun, to be our desired fitting
%function: 

ft = fittype('gfun(t0,tmax,ymax,alpha,time)', 'independent', {'time'},'coefficients', {'ymax', 'alpha','t0','tmax'}); 

%obtain fitting parameters: We have two different files, one human, one mouse 
%conditional below specifies this: 

if strcmp(species,'mouse')
prefs = parse_preference_file('mouse_AIF_fit_prefs.txt',1,...
    {'aif_lower_limits' 'aif_upper_limits' 'aif_initial_values' ...
     'aif_TolX' 'aif_MaxIter' 'aif_MaxFunEvals' 'aif_Robust' 'aif_TolFun'});
    prefs.aif_lower_limits = str2num(prefs.aif_lower_limits); 
    prefs.aif_upper_limits = str2num(prefs.aif_upper_limits);
    prefs.aif_initial_values = str2num(prefs.aif_initial_values); 
    prefs.aif_TolX = str2num(prefs.aif_TolX); 
    prefs.aif_MaxIter = str2num(prefs.aif_MaxIter);
    prefs.aif_MaxFunEvals = str2num(prefs.aif_MaxFunEvals); 
    prefs.aif_TolFun = str2num(prefs.aif_TolFun); 
elseif strcmp(species,'human')
    prefs = parse_preference_file('human_AIF_fit_prefs.txt',1,...
    {'aif_lower_limits' 'aif_upper_limits' 'aif_initial_values' ...
     'aif_TolX' 'aif_MaxIter' 'aif_MaxFunEvals' 'aif_Robust' 'aif_TolFun'});
    prefs.aif_lower_limits = str2num(prefs.aif_lower_limits); 
    prefs.aif_upper_limits = str2num(prefs.aif_upper_limits);
    prefs.aif_initial_values = str2num(prefs.aif_initial_values); 
    prefs.aif_TolX = str2num(prefs.aif_TolX); 
    prefs.aif_MaxIter = str2num(prefs.aif_MaxIter);
    prefs.aif_MaxFunEvals = str2num(prefs.aif_MaxFunEvals); 
    prefs.aif_TolFun = str2num(prefs.aif_TolFun); 
end

% collecting the fit paramters into an options structure: 
options = fitoptions('Method', 'NonlinearLeastSquares',...
    'MaxIter', prefs.aif_MaxIter,...
    'MaxFunEvals', prefs.aif_MaxFunEvals,...
    'TolFun', prefs.aif_TolFun,...
    'TolX', prefs.aif_TolX,...
    'Display', 'off',...
    'Lower', prefs.aif_lower_limits, ...
    'Upper', prefs.aif_upper_limits,...
    'StartPoint', prefs.aif_initial_values);

%Performing the fit: 
time = time_vect;

clear fit; 
%[fit_out] = fit(time,meanAIF,ft,'Lower', [3 1.5 0 0.01 ],'Upper', [ 100 10 0.05 0.1 ], 'StartPoint', [ 5 1.8 0.01 0.05]); 
[fit_out] = fit(time,meanAIF,ft,options);
% plot(fit_out); 
% hold on; 
% plot(time,meanAIF);
% hold off; 

% Now we pull out the coefficients from the fitted model and generate a
% generate the curve: 
ymax = fit_out.ymax; 
alpha= fit_out.alpha; 
t0 = fit_out.t0; 
tmax = fit_out.tmax; 
kappa = 0.04;   % a constant.  

%Plug them into the modeling equation: 
gt = gfun(t0,tmax,ymax,alpha,time);
Ct = zeros(numel(gt),1); 

for i = 1 : numel(gt) 
    Ct(i) = gt(i) + kappa * trapz(time(1:i), gt(1:i),1); 
end 

figure; 
time_vect_sec = time_vect * 60; %get time in seconds to plot
time_vect_sec = time_vect_sec + (bolus_time);
%base_time_vect_sec = base_time_vect * 60;
plot(time_vect_sec,Ct,'g'); 
hold on; 
plot(time_vect_sec,meanAIF, 'm');
title('AIF and Fitted AIF Over time');
xlabel('Time (s)');                                                                             
ylabel('Concentration (mM)');
legend('AIF', 'Fitted AIF');
hold off; 

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
function import_aif_CreateFcn(hObject, eventdata, handles)
% hObject    handle to import_aif (see GCBO)
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
set(handles.import_aif,'Value',0); 
set(handles.add_AIF_roi,'enable', 'off'); 
set(handles.AIF_path,'enable', 'off'); 

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'auto_cluster'
        handles.AIF_type =  0; 
        set(handles.add_AIF_roi,'enable', 'off'); 
        set(handles.AIF_path,'enable', 'off'); 
        set(handles.auto_cluster,'Value',1); 
        set(handles.user_aif,'Value',0);
        set(handles.import_aif,'Value',0); 
    case 'user_aif'
       handles.AIF_type = 1; 
       set(handles.add_AIF_roi,'enable', 'on'); 
       set(handles.AIF_path,'enable', 'on');
       set(handles.auto_cluster,'Value',0); 
       set(handles.user_aif,'Value',1);
       set(handles.import_aif,'Value',0); 
    case 'import_aif'
       handles.AIF_type = 2; 
       set(handles.add_AIF_roi,'enable', 'off'); 
       set(handles.AIF_path,'enable', 'off');
       set(handles.import_aif,'Value',1); 
       set(handles.auto_cluster,'Value',0); 
       set(handles.user_aif,'Value',0); 
end
guidata(hObject, handles);
