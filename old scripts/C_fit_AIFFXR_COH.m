%% C_fit_AIFFXR_COH:
%{
Fit to model with Shutter speed
Required:

smoothcurve.m % Smooth out the tumor ROI curve
FXRAIFhelper.m % Wrapper function for SS fitting
FXRAIF.m 

Thomas Ng
Caltech
Dec 2011

%}

clear all
%close all

disp('FXR')

%% Toggle options

% Set to 1 if you wantto smooth out the tumor curve with a fit. Here we use
% a smoothing function defined by Brix et. al, but you can alter it.
fitted = 1; 

%Time limit (in minutes), 0 if all to include in the dataset fitting
%Starting time set
nstarter = 0;
ntimeend = 12.5;

% Define the time when the injection started (in minutes)
ninject = 2.5;


%% 1. Load files from A_make_R1maps


[gogo,PathName1,FilterIndex] = uigetfile(['/data/studies/' '/*R1info.mat'],'Choose R1 file'); 

% Load the data files
directory = PathName1
rootname  = strrep(gogo, '.nii', '');
load(fullfile(PathName1, gogo));

starter = nstarter;
timeend= ntimeend;
inject  = ninject;

Cp = Cp(:);

%% 2. Time management
% Convert time to minutes
timeres = timeres/60;

%Generate time curve
timer   = [0:timeres:timeres*(numel(Cp)-1)];


if(starter)
    
    timecheck = abs(timer-starter);
    starter = find(timecheck == min(timecheck));
else
   % timeend = numel(timer);
    starter = 1;
end

%
if(timeend)
    
    timecheck = abs(timer-timeend);
    timeend = find(timecheck == min(timecheck));
else
    timeend = numel(timer);
end

%% 3. ROI fitting

% Prepare for lsqcurvefit
xdata{1}.timer = timer(starter:timeend);
timer          = xdata{1}.timer;
CpROI          = mean(Cp,2);
xdata{1}.Cp    = CpROI(starter:timeend);
CpROI          = xdata{1}.Cp;
R1t            = mean(R1tTOI,2);

xdata{1}.inject= inject;

% Change the function if you want alternate smoothing of the curve
if(fitted)
    xdata{1}.R1t   = smoothcurve(R1t(starter:timeend), xdata, 1);
   
else
    xdata{1}.R1t   = R1t(starter:timeend);
    
    xdata{1}.R1tTOI= R1tTOI(starter:timeend,:);
    
  
end
R10            = mean(1./(T1TUM));
xdata{1}.R10   = mean(R10);
xdata{1}.R1i   = mean(R10);
xdata{1}.fw    = fw;
xdata{1}.r1    = r1;


R1tTOI = R1tTOI(starter:timeend,:);

tic
ROIX = FXRAIFhelper(xdata, 1)
toc

Ctfitted = FXRAIF(ROIX, xdata);

b = figure; plot(timer,R1t(starter:timeend),'r.'), hold on, plot(timer, Ctfitted,'b'), plot(timer,xdata{1}.R1t,'g'),hold off,
title(['2 compartment model FXR: Ktrans: ' num2str(ROIX(1)) ' ve: ' num2str(ROIX(2)) ' tau_{i}: ' num2str(ROIX(3))] ), ylabel('R1t (1/s)'), xlabel('time (min)')
saveas(b, fullfile(PathName1, [rootname 'FXRfit.fig']));
%% 4. Save for voxel by voxel fitting, setup ready to go

xdata{1}.timer = timer;
xdata{1}.R1t    = R1tTOI;%(starter:timeend,:);
numvoxels      = size(Ct,2);
%res = [0.25 0.25 2];
save(fullfile(PathName1, [rootname 'AIF_FXR_ROI.mat']), 'res','currentimg','dynamname', 'lvind', 'tumind', 'ROIX', 'xdata', 'numvoxels');

