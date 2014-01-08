function results = B_AIF_fitting_func(results_a_path,start_time,end_time,start_injection,end_injection,fit_aif,average_aif,time_resolution)

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
%  average_aif        - average two AIFs together (bool)
%  time_resolution    - time resolution of dynamic scan (in minutes)
%
%
% The script loads the data arrays generated from A_make_R1maps_func(). 
% Then one can 
% 
% a) fit to a model.
% b) limit the time interval from which further data fitting is performed.
% c) There is an option to average the AIF from multiple datasets should you
% desire. This may be useful if you want to get better SNR or comparing
% multiday studies.
% 
% Requires:
% timeresmaker.m
% setupxdata2.m
% setupxdata3.m
% prepandsave.m
% AIFbiexpfithelp.m
% 
% Thomas Ng 
% Caltech
% December 2011
% Updated Dec 2013 Sam Barnes

%% Toggle options
% Replace curve from other dataset
replace = 0;
% If you have a threshold value for noise (testing).
threshold = 0;

%% DO NOT ALTER BELOW UNLESS YOU KNOW WHAT YOU ARE DOING

%% 1. Load the data array from previous script
load(results_a_path);

% Log input results
log_path = fullfile(PathName1, ['B_' rootname 'fitted_R1info.log']);
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
disp('User selected average two AIFs: ');
disp(average_aif);
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

%% In case of averaging.

Cp1 = Cp;

%% 3. Load the second one if averaged

% Averaged not functional
if(average_aif)
    [gogoav2,PathName1av2,FilterIndex] = uigetfile(['/data/studies/' '/*R1info.mat'],'Choose second average R1 file');
    
    % Load the data files
    PathName1av2
    
    load(fullfile(PathName1av2, gogoav2));
    
    load('timercut.mat');
    timeresmaker %@@@@@@@@@@@@@ Replace with function on inline code!
    
    Cp2 = Cp;
    
    %% reload the first one
    load(results_a_path);
    
    load('timercut.mat');
    timeresmaker %@@@@@@@@@@@@@ Replace with function on inline code!
    
    
    Cp = [Cp1 Cp2];
	% Setup the xdata for processing is done with setupxdata and setupxdata3.
	% These scripts sets ups the dataset to allow lsqcurvefit to work.
    setupxdata3 %@@@@@@@@@@@@@ Replace with function on inline code!
else
    %% Setup data to allow lsqcurvefit to work
	xdata{1}.timer = timer(start_time:end_time)';
	timer          = xdata{1}.timer;
	CpROI          = mean(Cp,2);
	CpROI          = CpROI(start_time:end_time);
	xdata{1}.Cp    = CpROI';

	Cp1 = Cp1(start_time:end_time,:);
	%Cp2 = Cp2(starter:timeend,:);

	% threshold (Remove noise manually)
	if(threshold)
		ind = find(CpROI > threshold);

		for j = 1:numel(ind)
			CpROI(ind(j)) = [];
			timer(ind(j)) = [];
			Cp1(ind(j), :)= [];
			Cp2(ind(j),:) = [];
		end

		xdata{1}.Cp    = CpROI';
		xdata{1}.timer = timer;
	end
end



%% 4. Define the step function
xdata{1}.step = [start_injection end_injection];

%% 5. Either Fit the curve to the model, keep the curve as raw data, or replace the curve with an AIF from another dataset.
if(~replace)
    if(fit_aif)
        [ROIX xAIF xdataAIF] = AIFbiexpfithelp(xdata, 1);
        
%         Cp = CpROI;
%         size(ROIX)
        Cp = ROIX;
    else
        ROIX = CpROI;
        Cp = CpROI;
    end
else
    [gogoav3,PathName1av3,FilterIndex] = uigetfile(['/data/studies/' '/*R1info*.mat'],'Choose fitted R1 file');
    
    load(fullfile(PathName1av3, gogoav3), 'Cp');
    ROIX = Cp;
end

%Ctfitted = FXLStep1AIFcfit(ROIX(1), ROIX(2), xdata{1}.Cp, xdata{1}.timer);

b = figure;
plot(timer,CpROI,'r.');
hold on;
plot(timer, ROIX,'b');

M{1} = 'Plasma curve';
M{2} = 'Fitted curve';
if(average_aif)
    plot(timer, mean(Cp1,2), 'gx'), plot(timer, mean(Cp2,2), 'kx'),
    M{end+1} = 'Averaged curve 1';
    M{end+1} = 'Averaged curve 2';
end
legend(M);
root_modified = rootname;
root_modified(end)='';
hold off;
title([root_modified ' - AIF Bi-Exponential, Linear Upslope'], 'Interpreter', 'none');
ylabel('C_t (mM)');
xlabel('time (min)');
saveas(b, fullfile(PathName1, [rootname 'AIF_fitting.fig']));

%% 6.  Save for voxel by voxel fitting, setup ready to go
xdata{1}.timer = timer;
xdata{1}.Ct    = Ct(start_time:end_time,:);
numvoxels      = size(Ct,2);

save(fullfile(PathName1, ['B_' rootname 'fitted_R1info.mat']));

if(average_aif)
    % Save the fitted
    save('temp.mat', 'Cp', 'PathName1av2', 'gogoav2', 'threshold', 'Choose', 'Averaged', 'fit', 'edit');
    
    clear all
    
    load('temp.mat');
    % Load the data files
    directory = PathName1av2;
    rootname  = strrep(gogoav2, '.nii', '');
    load(fullfile(PathName1av2, gogoav2));
    load('timercut.mat');
    timeresmaker
    setupxdata2
    load('temp.mat');
    prepandsave
    disp('Processed 2nd timepoint');
end

results = fullfile(PathName1, ['B_' rootname 'fitted_R1info.mat']);
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
