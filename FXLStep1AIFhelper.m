%% Helper Function for nonlinear curvefit to FXLAIF model, no vp
function x = FXLStep1AIFhelper(Ct,Cp,timer)
  
warning off
% Ct = xdata{1}.Ct;
% i  = voxel;
% Ct = Ct(:,i);

% tic

% %configure the optimset for use with lsqcurvefit
% options = optimset('lsqcurvefit');
% 
% %increase the number of function evaluations for more accuracy
% options.MaxFunEvals = 50;
% options.MaxIter     = 50;
% options.TolFun      = 10^(-12);
% options.TolX        = 10^(-6);
% options.Diagnostics = 'off';
% options.Display     = 'off';
% options.Algorithm   = 'levenberg-marquardt';
% 
% % lb = [0 0.03];
% % ub = [Inf 1];
% lb = [0 10^20];
% ub = [Inf Inf];
% 
% [x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@FXLStep1AIF, ...
%     [0.005 0.01], xdata, ...
%     Ct',lb,ub,options);
% 
% %x(1) =			% ktrans
% %x(2) = 		% ve
% x(3) = 1;		% unused
% x(4) = resnorm;	% residual

% Get values from pref file
prefs = parse_preference_file('dce_preferences.txt',0,...
	{'voxel_lower_limit_ktrans' 'voxel_upper_limit_ktrans' 'voxel_initial_value_ktrans' ...
	'voxel_lower_limit_ve' 'voxel_upper_limit_ve' 'voxel_initial_value_ve' ...
	'voxel_TolFun' 'voxel_TolX' 'voxel_MaxIter' 'voxel_MaxFunEvals' 'voxel_Robust'});
lower_limit_ktrans = str2num(prefs.voxel_lower_limit_ktrans);
upper_limit_ktrans = str2num(prefs.voxel_upper_limit_ktrans);
initial_value_ktrans = str2num(prefs.voxel_initial_value_ktrans);
lower_limit_ve = str2num(prefs.voxel_lower_limit_ve);
upper_limit_ve = str2num(prefs.voxel_upper_limit_ve);
initial_value_ve = str2num(prefs.voxel_initial_value_ve);
TolFun = str2num(prefs.voxel_TolFun);
TolX = str2num(prefs.voxel_TolX);
MaxIter = str2num(prefs.voxel_MaxIter);
MaxFunEvals = str2num(prefs.voxel_MaxFunEvals);
Robust = prefs.voxel_Robust;

% Use Curvefitting tool box instead of optimization toolbox (lsqcurvefit)
% as curvefitting will easily return confidence intervals on the fit
% performance of the two appears to be the same
options = fitoptions('Method', 'NonlinearLeastSquares',...
	'Algorithm', 'Levenberg-Marquardt',...
	'MaxIter', MaxIter,...
	'MaxFunEvals', MaxFunEvals,...
	'TolFun', TolFun,...
	'TolX', TolX,...
	'Display', 'off',...
	'Lower',[lower_limit_ktrans lower_limit_ve],...
	'Upper', [upper_limit_ktrans upper_limit_ve],...
	'StartPoint', [initial_value_ktrans initial_value_ve],...
	'Robust', Robust);
ft = fittype('FXLStep1AIFcfit( Ktrans, ve, Cp, T1)',...
	'independent', {'T1', 'Cp'},...
	'coefficients',{'Ktrans', 've'});
[f, gof, output] = fit([timer, Cp'],Ct,ft, options);
confidence_interval = confint(f,0.95);
% toc

%Calculate the R2 fit
x(1) = f.Ktrans;			% ktrans
x(2) = f.ve;				% ve
x(3) = 1;					% unused
x(4) = gof.sse;				% residual
% x(5) = confidence_interval;	% 95 CI
x(5) = confidence_interval(1,1);% (95 lower CI of ktrans)
x(6) = confidence_interval(2,1);% (95 upper CI of ktrans)
x(7) = confidence_interval(1,2);% (95 lower CI of ve)
x(8) = confidence_interval(2,2);% (95 upper CI of ve)
