%% This file is used to apply a model fitting to the arterial input function.
%{

The script loads the data arrays generated from the previous step. 
Then, if desired, one can 

a) fit to a model.
b) limit the time interval from which further data fitting is performed.
c) There is an option to average the AIF from multiple datasets should you
desire. This may be useful if you want to get better SNR or comparing
multiday studies.

Requires:
timeresmaker.m
setupxdata.m
setupxdata3.m
prepandsave.m
AIFbiexpfithelp.m

Thomas Ng 
Caltech
December 2011

%}

function results = B_AIF_fitting_func(results_a_path,start_time,end_time,start_injection,end_injection,fit_aif,average_aif)

disp('Fitting AIF')

%% Toggle options
% 
% %Average between 2 curves (Will prompt you to choose the second dataarray
% %mat file)
% average_aif= 0;
% 
% % Do you want to do a model fit? (Set to 0 if you are doing reference
% % region)
% fit_aif = 1;
% 
% % The next two toggles allow you to choose a sub-interval from which to
% % fit.
% %Starting time set (in minutes, 0 if all)
% start_time = 0;
% %Time limit (in minutes), 0 if all
% end_time = 0;
% 
% % The next two toggles define the start and end of the injection duration
% % (for convolution fitting) in minutes.
% start_injection = 5.0;
% end_injection = 5.25; 

% Replace curve from other dataset
replace = 0;
% If you have a threshold value for noise (testing).
threshold = 0;

%% DO NOT ALTER BELOW UNLESS YOU KNOW WHAT YOU ARE DOING

%% 1. Load the data array from previous script

%     [gogoav1,PathName1av1,FilterIndex] = uigetfile(['C:\Users\sbarnes\Documents\data\6 DCE Stroke\sb01_06nov13.mH1' '/*R1info.mat'],'Choose original R1 file');
%     
%     
%     % Load the data files
%     directory = PathName1av1
%     rootname  = strrep(gogoav1, '.nii', '');
%     load(fullfile(PathName1av1, gogoav1));

[PathName1av1,base,ext] = fileparts(results_a_path);
gogoav1 = [base ext];
load(results_a_path);


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
    load(fullfile(PathName1av1, gogoav1));
    
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
    tic
    if(fit_aif)
        [ROIX xAIF xdataAIF] = AIFbiexpfithelp(xdata, 1);
        
        %Cp = CpROI;
        size(ROIX)
        Cp = ROIX;
    else
        ROIX = CpROI;
        Cp = CpROI;
    end
    toc
else
    [gogoav3,PathName1av3,FilterIndex] = uigetfile(['/data/studies/' '/*R1info*.mat'],'Choose fitted R1 file');
    
    load(fullfile(PathName1av3, gogoav3), 'Cp');
    ROIX = Cp;
end

%Ctfitted = FXLStep1AIFcfit(ROIX(1), ROIX(2), xdata{1}.Cp, xdata{1}.timer);

b = figure;plot(timer,CpROI,'r.'), hold on, plot(timer, ROIX,'b'),

M{1} = 'Plasma curve';
M{2} = 'Fitted curve';
if(average_aif)
    plot(timer, mean(Cp1,2), 'gx'), plot(timer, mean(Cp2,2), 'kx'),
    M{end+1} = 'Averaged curve 1';
    M{end+1} = 'Averaged curve 2';
    
end
legend(M)
hold off,
title([rootname ' - AIF Bi-Exponential, Linear Upslope'], 'Interpreter', 'none'), ylabel('C_t (mM)'), xlabel('time (min)')
saveas(b, fullfile(PathName1av1, [rootname 'AIF_fitting.fig']));

%% 6.  Save for voxel by voxel fitting, setup ready to go
xdata{1}.timer = timer;
xdata{1}.Ct    = Ct(start_time:end_time,:);
numvoxels      = size(Ct,2);

save(fullfile(PathName1, [rootname 'fitted_R1info.mat']));

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


results = fullfile(PathName1, [rootname 'fitted_R1info.mat']);

disp('Finished B');
