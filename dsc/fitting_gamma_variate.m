function [ Ct ] = fitting_gamma_variate( meanAIF, time_vect )
%FITTING_GAMMA_VARIATE Now we fit the AIF with a SCR model: 

%assigning the gamma variate function, gfun, to be our desired fitting
%function: 

ft = fittype('gfun(t0,tmax,ymax,alpha,time)', 'independent', {'time'},'coefficients', {'ymax', 'alpha','t0','tmax'}); 

%obtain fitting parameters:
prefs = parse_preference_file('dsc_preferences.txt',1,...
    {'aif_gamma_lower_limits' 'aif_gamma_upper_limits' 'aif_gamma_initial_values' ...
     'aif_TolX' 'aif_MaxIter' 'aif_MaxFunEvals' 'aif_Robust' 'aif_TolFun'});
    prefs.aif_gamma_lower_limits = str2num(prefs.aif_gamma_lower_limits); 
    prefs.aif_gamma_upper_limits = str2num(prefs.aif_gamma_upper_limits);
    prefs.aif_gamma_initial_values = str2num(prefs.aif_gamma_initial_values); 
    prefs.aif_TolX = str2num(prefs.aif_TolX); 
    prefs.aif_MaxIter = str2num(prefs.aif_MaxIter);
    prefs.aif_MaxFunEvals = str2num(prefs.aif_MaxFunEvals); 
    prefs.aif_TolFun = str2num(prefs.aif_TolFun); 

% collecting the fit paramters into an options structure: 
options = fitoptions('Method', 'NonlinearLeastSquares',...
    'MaxIter', prefs.aif_MaxIter,...
    'MaxFunEvals', prefs.aif_MaxFunEvals,...
    'TolFun', prefs.aif_TolFun,...
    'TolX', prefs.aif_TolX,...
    'Display', 'off',...
    'Lower', prefs.aif_gamma_lower_limits, ...
    'Upper', prefs.aif_gamma_upper_limits,...
    'StartPoint', prefs.aif_gamma_initial_values);

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

end

