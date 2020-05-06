% AIFbixexpfithelper is a wrapper function to fit the AIF to a biexponential model

%{

AIFbiexpcon.m is the file defining the function. You can of course alter
this as you wish

out is the fitted timecurve over the desired time interval
x is the parameters of the fit
xdata stores the input parameters of the fit (cp (meanAIF?? but over the whole run and timer can be passed in separately) 

the text file is parseed 
out is fitted concentration curve 
x data structure that holds the fitted unknowns 
verbose should be able to cut out as it is for testing purposes
%}

function [Cp, x, xdata, rsquare] = Forced_linear_AIFbiexpfithelplocal(xdata, verbose)

warning off

%{
if ~iscell(xdata)
    foo{1} = xdata;
    xdata = foo;
end
%}

Cp = xdata.Cp;
Cp = Cp(:);

t  = xdata.timer;
oldt= t;
oldt = oldt(:);

xdata.timer = t;
t  = t(:);

% Get preferences
prefs = parse_preference_file('dsc_preferences.txt',0,...
    {'aif_lower_limits' 'aif_upper_limits' 'aif_initial_values' ...
    'aif_TolFun' 'aif_TolX' 'aif_MaxIter' 'aif_MaxFunEvals' 'aif_Robust'});
lower_limits = str2num(prefs.aif_lower_limits);
upper_limits = str2num(prefs.aif_upper_limits);
initial_values = str2num(prefs.aif_initial_values);
TolFun = str2num(prefs.aif_TolFun);
TolX = str2num(prefs.aif_TolX);
MaxIter = str2num(prefs.aif_MaxIter);
MaxFunEvals = str2num(prefs.aif_MaxFunEvals);
Robust = prefs.aif_Robust;

if verbose>0
    fprintf('lower_limits = %s\n',num2str(lower_limits));
    fprintf('upper_limits = %s\n',num2str(upper_limits));
    fprintf('initial_values = %s\n',num2str(initial_values));
    fprintf('TolFun = %s\n',num2str(TolFun));
    fprintf('TolX = %s\n',num2str(TolX));
    fprintf('MaxIter = %s\n',num2str(MaxIter));
    fprintf('MaxFunEvals = %s\n',num2str(MaxFunEvals));
    fprintf('Robust = %s\n\n',Robust);
end

%configure the optimset for use with lsqcurvefit
options = optimset('lsqcurvefit');

%increase the number of function evaluations for more accuracy
options.MaxFunEvals = MaxFunEvals;
options.MaxIter     = MaxIter;
options.TolFun      = TolFun;
options.TolX        = TolX;
options.Diagnostics = 'off';
options.Display     = 'off';
options.Algorithm   = 'levenberg-marquardt';
options.Robust      = Robust;

% Choose upper and lower bounds only for trust-region methods.
% lb = [0 0 0 0];
% ub = [5 5 5 5];
% initial_values = [1 1 1 1];

%% Split the fitting between the biexponential phase and the linear phase
t = oldt;
if verbose>0
    figure, plot(t, Cp./max(Cp), 'b.');
    title('Weighting, Injection Period, AIF Curve'), xlabel('time (min)');
end


%[x y] = ginput(1);
% x = 1;
% temp = abs(x-t);
% ind  = find(temp == min(temp));

timer = t;
start = xdata.step;
ended = start(2);
start = start(1);

%{
start_index = find(abs(timer - start) == min(abs(timer - start)));
end_index   = find(abs(timer - ended) == min(abs(timer - ended)));
%}

step = zeros(size(timer));
step(start:ended) = 1; %originally step(start_index:end_index) =1;
xdata.step = step;


%W is the weighting matrix, should you want to emphasise certain
% datapoints
W = ones(size(Cp));

%calculate the most likely first peak by finding the first local maxima
%that is within 3% of the max value of the whole data set
[local_maxima, maxima_indexes] = findpeaks(Cp);
maxima_iterator = 1;
while(local_maxima(maxima_iterator) < (0.97 * max(local_maxima)))
    maxima_iterator = maxima_iterator + 1;
end

max_index = maxima_indexes(maxima_iterator);

% WW= sort(Cp.*step, 'descend');
% ind(1) = find(Cp == WW(1));
% ind(2) = find(Cp == WW(2));
% ind(3) = find(Cp == WW(3));

step(max_index+1:end) = 0;
xdata.step = step;

%max_index = ended; %this is a new addition which has made the function run
%in the past

if isempty(find(step==1,1))
    % Something has gone wrong, reset to default
    step = zeros(size(timer));
    step(start_index:end_index) = 1;
    xdata.step = step;
end

if verbose>0
    hold on, 
    plot(t,step, 'r'),
    plot(t(max_index), Cp(max_index)/max(Cp), 'kx', 'MarkerSize', 30);
end

% Alter the weightings here.
% W(max_index) =1;
% W(max_index+1)= 1;
% W(max_index-1)= 1;

if verbose>0
    plot(t, W, 'gx');
end

maxer = Cp(max_index);
xdata.maxer = maxer;
Cp = Cp.*W;



%end_baseline = find(xdata.step > 0);
%baseline = mean(Cp(1:end_baseline(1)));
%xdata{1}.baseline = baseline;
upper_limits(1) = maxer*2;
upper_limits(2) = maxer*2;
initial_values(1) = maxer*0.5;
initial_values(2) = maxer*0.5;
initial_values(3) = 1; 
initial_values(4) = 1; 
if verbose>0
    disp('Fitting raw values, limits and initial values adjusted');
    fprintf('lower_limits = %s\n',num2str(lower_limits));
    fprintf('upper_limits = %s\n',num2str(upper_limits));
    fprintf('initial_values = %s\n',num2str(initial_values));
end


% Currently, we use AIF
[x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@Forced_linear_AIFbiexpconlocal, ...
    initial_values, xdata, ...
    Cp',lower_limits,upper_limits,options);

xdata.timer = oldt;
rsquare = 1 - resnorm / norm(Cp-mean(Cp))^2;
if verbose>0
    disp(['R^2 of AIF fit = ' num2str(rsquare)]);
end

Cp = Forced_linear_AIFbiexpconlocal(x, xdata);



