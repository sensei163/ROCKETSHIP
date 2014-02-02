function results = B_AIF_fitting_func(results_a_path,start_time,end_time,start_injection,end_injection,fit_aif,import_aif_path,time_resolution)

% B_AIF_fitting_func - This file is used to apply a model fitting to the
% arterial input function. AIF is fit to a bi-exponential model with a
% linear ramp up (from baseline to peak)
% 
% Inputs:
%  results_a_path     - *.mat Results from part A
%  start_time         - start time to restrict analysis to (in minutes, 0 
%                       if all)
%  end_time           - end time to restrict analysis to (in minutes, 0 if 
%                       all)
%  start_injection    - start of contrast injection (in minutes)
%  end_injection      - end of contrast injection (in minutes)
%  fit_aif            - fit AIF to model (bool, 0 for reference region)
%  import_aif_path    - import AIF from .mat file (string)
%  time_resolution    - time resolution of dynamic scan (in minutes)
%
%
% The script loads the data arrays generated from A_make_R1maps_func(). 
% Then one can 
% 
% a) fit to a model. 
% b) limit the time interval from which further data fitting is performed. 
% c) There is an option import a population or averaged AIF from multiple
% datasets should you desire. This may be useful if you want to get better
% SNR or comparing multiday studies.
% 
% Requires:
% B_AIF_fitting_func.m
% AIFbiexpfithelp.m
% DataHash.m
% AIFbiexpcon.m
% parse_preference_file.m
% timeresmaker.m *
% setupxdata2.m *
% setupxdata3.m *
% prepandsave.m *
% 
% Thomas Ng 
% Caltech
% December 2011
% Updated Dec 2013 Sam Barnes

%% Toggle options
% If you have a threshold value for noise (testing).
threshold = 0;

%% DO NOT ALTER BELOW UNLESS YOU KNOW WHAT YOU ARE DOING

%% 1. Load the data array from previous script

load(results_a_path);

% unload the variables from previous data array

rootname = Adata.rootname;
Cp       = Adata.Cp;
Ct       = Adata.Ct;


% update output path to be same as location of input
[PathName1,~,~] = fileparts(results_a_path);

% Log input results
log_path = fullfile(PathName1, ['B_' rootname 'test_R1info.log']);
if exist(log_path, 'file')==2
  delete(log_path);
end
diary(log_path);
fprintf('************** User Input **************\n\n');
disp('User selected part A results: ');
disp(results_a_path);
Opt.Input = 'file';
a_md5 = DataHash(results_a_path, Opt);
fprintf('File MD5 hash: %s\n\n', a_md5)
disp('User selected start time (min): ');
disp(start_time);
disp('User selected end time (min): ');
disp(end_time);
disp('User selected start injection (min): ');
disp(start_injection);
disp('User selected end injection (min): ');
disp(end_injection);
disp('User selected fit AIF to model: ');
disp(fit_aif);
disp('User selected import AIF path: ');
fprintf('%s\n\n',import_aif_path);
disp('User selected time resolution (sec)');
disp(time_resolution*60);
fprintf('************** End User Input **************\n\n\n');

disp('Starting Part B - Fitting AIF')
disp(datestr(now))
disp(' ');
tic

%% 2. time management (Change the starter, timeend if needed)
%% b) Time management

%Generate time curve
timer   = 0:time_resolution:time_resolution*(size(Cp,1)-1);

%Restrict to specified times
if(start_time)
    timecheck = abs(timer-start_time);
    start_time = find(timecheck == min(timecheck));
else
    start_time = 1;
end

if(end_time)
    timecheck = abs(timer-end_time);
    end_time = find(timecheck == min(timecheck));
else
    end_time = numel(timer);
end


%% Setup data to allow lsqcurvefit to work
xdata{1}.timer = timer(start_time:end_time)';
timer          = xdata{1}.timer;
CpROI          = mean(Cp,2);
CpROI          = CpROI(start_time:end_time)';
xdata{1}.Cp    = CpROI;

% threshold (Remove noise manually)
if(threshold)
    ind = find(CpROI > threshold);

    for j = 1:numel(ind)
        CpROI(ind(j)) = [];
        timer(ind(j)) = [];
        Cp1(ind(j), :)= [];
        Cp2(ind(j),:) = [];
    end

    xdata{1}.Cp    = CpROI;
    xdata{1}.timer = timer;
end




%% 4. Define the step function
xdata{1}.step = [start_injection end_injection];

%% 5. Either Fit the curve to the model, keep the curve as raw data, or replace the curve with an AIF from another dataset.
M{1} = '';
aif_name = '';
if isempty(import_aif_path)
    if(fit_aif)
        [Cp_fitted xAIF xdataAIF] = AIFbiexpfithelp(xdata, 1);
        Cp_use = Cp_fitted;
        M{2} = 'Fitted Curve';
        aif_name = 'fitted';
    else
        Cp_use = CpROI;
        M{2} = 'Using Raw Curve';
        aif_name = 'raw';
    end
else
    external = load(import_aif_path);
    if isfield(external,'Cp_use')
        Cp_use = external.Cp_use;
    elseif isfield(external,'xdata')
        Cp_use = external.xdata{1}.Cp;
    else
        disp('No Cp curve found in selected file')
        return
    end
    M{2} = 'Imported Curve';
    aif_name = 'imported';
end

% 5.5 Plot the results
b = figure;
plot(timer,CpROI,'r.');
hold on;
plot(timer, Cp_use,'b');

M{1} = 'Original Plasma Curve';
% M{2} = 'Selected Curve';

legend(M);
root_modified = rootname;
root_modified(end)='';
hold off;
title([root_modified ' - AIF Bi-Exponential, Linear Upslope'], 'Interpreter', 'none');
ylabel('Concentration (mM)');
xlabel('Time (min)');
saveas(b, fullfile(PathName1, [rootname 'AIF_fitting.fig']));

%% 6.  Save for voxel by voxel fitting, setup ready to go
xdata{1}.timer = timer;
xdata{1}.Ct    = Ct(start_time:end_time,:);
numvoxels      = size(Ct,2);
xdata{1}.Cp = Cp_use;

%% 7. Setup  Bresults to output
Bdata.CpROI         = CpROI;
Bdata.Cp_use        = Cp_use;
Bdata.aif_name      = aif_name;
Bdata.end_injection = end_injection;
Bdata.end_time      = fit_aif;
Bdata.log_path      = log_path;
Bdata.numvoxels     = numvoxels;
Bdata.results_a_path= results_a_path; % Location of associated A data
Bdata.start_injection=start_injection;
Bdata.start_time    = start_time;
Bdata.threshold     = threshold;
Bdata.time_resolution=time_resolution;
Bdata.timer         = timer;
Bdata.xdata         = xdata;

results = fullfile(PathName1, ['B_' rootname aif_name '_R1info.mat']);

save(results, 'Bdata');

Opt.Input = 'file';
mat_md5 = DataHash(results, Opt);
disp(' ')
disp('MAT results saved to: ')
disp(results)
disp(['File MD5 hash: ' mat_md5])

disp(' ');
disp('Finished B');
disp(datestr(now))
toc
diary off;

% Keep log consistent naming with data file
movefile(fullfile(PathName1, ['B_' rootname 'test_R1info.log']) , fullfile(PathName1, ['B_' rootname aif_name '_R1info.log']));
