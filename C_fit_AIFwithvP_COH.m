%% C_fit_AIFFXL_COH:
%{
Fit to model with vp SM2
Required:

smoothcurve.m % Smooth out the tumor ROI curve
FXLStep1AIFhelper_vp.m % Wrapper function for SM2 fitting

Thomas Ng
Caltech
Dec 2011

%}

clear all
%close all
disp('FXL')

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
    %timeend = numel(timer);
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
CtROI          = mean(Ct,2);
oCtROI         = CtROI(starter:timeend);
xdata{1}.inject= inject;

% Change the function if you want alternate smoothing of the curve
if(fitted)
    xdata{1}.Ct= smoothcurve(CtROI(starter:timeend), xdata, 0);
else

xdata{1}.Ct    = CtROI(starter:timeend);
end
CtROI          = xdata{1}.Ct;

tic
ROIX = FXLStep1AIFhelper_vp(xdata,1)
toc


%% Multisearch: This may be useful for global optimization. Commented out for now, but uncomment to use.
% matlabpool local 7
% RS = RandomStartPointSet('NumStartPoints', 50);
% options = optimset('Algorithm','trust-region-reflective');
% problem = createOptimProblem('lsqcurvefit', 'x0', ROIX(1:3), 'objective',@FXLStep1AIF_vp, 'lb', [0 0 0], 'ub', [1 1 1], 'xdata', xdata, 'options', options, 'ydata', CtROI' );
% ms = MultiStart('StartPointsToRun', 'bounds', 'UseParallel', 'always', 'Display', 'off');
% [xming,fming,flagg,outptg,manyminsg] = run(ms,problem, RS);
% disp('Applying Multistart Global optimization KtransRR, Ktrans TOI, ve TOI:'),
% ROIX = xming
% matlabpool close
%%

Ctfitted = FXLStep1AIF_vp(ROIX,xdata);

b = figure; plot(timer,oCtROI,'r.'), hold on, plot(timer, Ctfitted,'b'),plot(timer,CtROI,'g'), hold off, 
title(['2 compartment model vp: Ktrans: ' num2str(ROIX(1)) ' ve: ' num2str(ROIX(2)) ' vp: ' num2str(ROIX(3))] ), ylabel('Ct (mM)'), xlabel('time (min)')
saveas(b, fullfile(PathName1, [rootname 'FXL_vp_fit.fig']));
%% 4. Save for voxel by voxel fitting, setup ready to go

xdata{1}.timer = timer;
xdata{1}.Ct    = Ct(starter:timeend,:);
numvoxels      = size(Ct,2);
%res = [0 0.25 0.25 2];
save(fullfile(PathName1, [rootname 'AIF_with_vpFIT_ROI.mat']), 'res', 'currentimg', 'dynamname', 'lvind', 'tumind', 'ROIX', 'xdata', 'numvoxels');
        