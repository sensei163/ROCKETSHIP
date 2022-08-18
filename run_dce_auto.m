function run_dce_auto(subject_tp_path)
% Find and add subpaths 
 mfilepath=fileparts(which('run_dce_auto'));
 addpath(fullfile(mfilepath,'dce'));
 addpath(fullfile(mfilepath,'external_programs'));
 addpath(fullfile(mfilepath,'external_programs/niftitools'));
 addpath(fullfile(mfilepath,'parametric_scripts'));
%% RUN A
% load A prefs
script_prefs = parse_preference_file('script_preferences.txt', 0, ...
    {'noise_pathpick' 'noise_pixsize' 'dynamic_files' ...
    'aif_files' 'roi_files' 't1map_files' 'noise_files' 'drift_files' ...
    'rootname' 'fileorder' 'quant' 'roimaskroi' 'aifmaskroi' 'aif_rr_type' ...
    'tr' 'fa' 'hematocrit' 'snr_filter' 'relaxivity' 'injection_time' ...
    'injection_duration' 'drift_global' 'blood_t1'});

% force 4D files
filevolume = 1;
% don't need to display pretty file list
LUT = 1;

% type casts
noise_pathpick = str2num(script_prefs.noise_pathpick);
noise_pixsize = str2num(script_prefs.noise_pixsize);

dynamic_files = cellstr(strcat(subject_tp_path,script_prefs.dynamic_files));
aif_files = cellstr(strcat(subject_tp_path,script_prefs.aif_files));
roi_files = cellstr(strcat(subject_tp_path,script_prefs.roi_files));
t1map_files = cellstr(strcat(subject_tp_path,script_prefs.t1map_files));

if ~strcmp(script_prefs.noise_files,'')
    noise_files = cellstr(strcat(subject_tp_path,script_prefs.noise_files));
else
    noise_files = '';
end;

if ~strcmp(script_prefs.drift_files,'')
    drift_files = cellstr(strcat(subject_tp_path,script_prefs.drift_files));
else
    drift_files = '';
end;
    
quant = str2num(script_prefs.quant);
roimaskroi = str2num(script_prefs.roimaskroi);
aifmaskroi = str2num(script_prefs.aifmaskroi);
tr = str2double(script_prefs.tr);
fa = str2double(script_prefs.fa);
hematocrit = str2double(script_prefs.hematocrit);
snr_filter = str2num(script_prefs.snr_filter);
injection_time = str2num(script_prefs.injection_time);
relaxivity = str2double(script_prefs.relaxivity);
drift_global = str2num(script_prefs.drift_global);
blood_t1 = str2num(script_prefs.blood_t1);
injection_duration = str2num(script_prefs.injection_duration);

% main function call
[A_results, errormsg] = A_make_R1maps_func(filevolume, noise_pathpick, ...
    noise_pixsize, LUT, dynamic_files,aif_files, ...
    roi_files, t1map_files, ...
    noise_files, drift_files, ...
    script_prefs.rootname, script_prefs.fileorder, quant, roimaskroi, ...
    aifmaskroi, script_prefs.aif_rr_type, tr, fa, hematocrit, snr_filter, ...
    relaxivity, injection_time, drift_global, blood_t1, injection_duration);

if ~isempty(errormsg)
    
    disp_error(errormsg, handles);
    return;
end

%% RUNB
% load B prefs
script_prefs = parse_preference_file('script_preferences.txt', 0, ...
    {'start_time', 'end_time', 'auto_find_injection', 'start_injection', ...
    'end_injection', 'fit_aif', 'time_resolution', 'aif_type', ...
    'import_aif_path', 'timevectyn', 'timevectpath'});

% type casts
start_time = str2num(script_prefs.start_time);
end_time = str2num(script_prefs.end_time);
auto_find_injection = str2num(script_prefs.auto_find_injection);
start_injection = str2num(script_prefs.start_injection);
end_injection = str2num(script_prefs.end_injection);
%fit_aif = str2num(script_prefs.fit_aif);
time_resolution = str2double(script_prefs.time_resolution);
aif_type = str2num(script_prefs.aif_type);
import_aif_path = script_prefs.import_aif_path;
timevectyn = str2num(script_prefs.timevectyn);
timevectpath = script_prefs.timevectpath;

% convert time resolution into minutes
time_resolution = time_resolution / 60;

fit_aif = aif_type;
if (auto_find_injection)
    start_injection = -1;
    end_injection = -1;
end

% main function call
B_results = B_AIF_fitting_func(A_results,start_time,end_time, start_injection,end_injection,fit_aif,import_aif_path,time_resolution, timevectpath);

%% RUND
% whatever happened to C?

