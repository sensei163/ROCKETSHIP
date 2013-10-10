%% This file is used to apply a model fitting to the arterial input function.
disp('Fitting AIF')
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

%% Clear workspace
clear all
%close all

%% Toggle options

%Average between 2 curves (Will prompt you to choose the second dataarray
%mat file)
Averaged= 0;

% Do you want to do a model fit? (Set to 0 if you are doing reference
% region)
fit = 1;

% The next two toggles allow you to choose a sub-interval from which to
% fit.
%Starting time set (in minutes, 0 if all)
starter = 0;
%Time limit (in minutes), 0 if all
timeend = 25;

% The next two toggles define the start and end of the injection duration
% (for convolution fitting) in minutes.
start = 2.5;
ended = 2.7; 

% For testing only
edit= 0;

% If you have a threshold value for noise (testing).
threshold = 0;
% Replace curve from other dataset
replace = 0;

%% DO NOT ALTER BELOW UNLESS YOU KNOW WHAT YOU ARE DOING

%% 1. Load the data array from previous script

    [gogoav1,PathName1av1,FilterIndex] = uigetfile(['/data/studies/' '/*R1info.mat'],'Choose original R1 file');
    
    
    % Load the data files
    directory = PathName1av1
    rootname  = strrep(gogoav1, '.nii', '');
    load(fullfile(PathName1av1, gogoav1));

% if(size(Cp,2) < size(Cp,1))
%   Cp = Cp';
% end


%% 2. time management (Change the starter, timeend if needed)

timeresmaker

%% In case of averaging.

Cp1 = Cp;

%% 3. Load the second one if averaged

if(Averaged)
    [gogoav2,PathName1av2,FilterIndex] = uigetfile(['/data/studies/' '/*R1info.mat'],'Choose second average R1 file');
    
    
    % Load the data files
    PathName1av2
    
    load(fullfile(PathName1av2, gogoav2));
    
    load('timercut.mat');
    timeresmaker
    
    Cp2 = Cp;
    
    %% reload the first one
    load(fullfile(PathName1av1, gogoav1));
    
    load('timercut.mat');
    timeresmaker
    
    
    Cp = [Cp1 Cp2];
    setupxdata3
else
    setupxdata
end


% Setup the xdata for processing is done with setupxdata and setupxdata3.
% These scripts sets ups the dataset to allow lsqcurvefit to work.


%% 4. Define the step function



% starter1 = find((timer - start) == min(abs(timer - start)));
% ender   = find((timer - ended) == min(abs(timer - ended)));
%
% step = zeros(size(timer));
% step(starter1:ender) = 1;

xdata{1}.step = [start ended];

%% 5. Either Fit the curve to the model, keep the curve as raw data, or replace the curve with an AIF from another dataset.
if(~replace)
    tic
    if(fit)
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
if(Averaged)
    plot(timer, mean(Cp1,2), 'gx'), plot(timer, mean(Cp2,2), 'kx'),
    M{end+1} = 'Averaged curve 1';
    M{end+1} = 'Averaged curve 2';
    
end
legend(M)
hold off,
title([rootname 'AIF fitting with bi exponential, linear upslope'] ), ylabel('Ct (mM)'), xlabel('time (min)')
saveas(b, fullfile(PathName1av1, [rootname 'AIF_fitting.fig']));

%% 6.  Save for voxel by voxel fitting, setup ready to go
prepandsave

if(Averaged)
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


filename = fullfile(PathName1, [rootname '_fitted_R1info.mat']);


