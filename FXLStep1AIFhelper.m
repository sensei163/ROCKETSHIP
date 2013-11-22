%% Helper Function for nonlinear curvefit to FXLAIF model, no vp
function x = FXLStep1AIFhelper(xdata, voxel,smooth_model,smoothing_window)
if nargin < 3
    smooth_model = 0;
end
  
warning off
Ct = xdata{1}.Ct;
i  = voxel;
Ct = Ct(:,i);

t = xdata{1}.timer;
xdata{1}.timer = t(:);

if strcmp(smooth_model,'moving')
	Ct = smooth(Ct,smoothing_window,'moving');
elseif strcmp(smooth_model,'rlowess')
	Ct = smooth(Ct,smoothing_window/size(Ct,1),'rlowess');
else
	% no smoothing
end

%configure the optimset for use with lsqcurvefit
options = optimset('lsqcurvefit');

%increase the number of function evaluations for more accuracy
options.MaxFunEvals = 50;
options.MaxIter     = 50;
options.TolFun      = 10^(-12);
options.TolX        = 10^(-6);
options.Diagnostics = 'off';
options.Display     = 'off';
options.Algorithm   = 'levenberg-marquardt';

% lb = [0 0.03];
% ub = [Inf 1];
lb = [0 10^20];
ub = [Inf Inf];

% tic
[x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@FXLStep1AIF, ...
    [0.005 0.01], xdata, ...
    Ct',lb,ub,options);

%x(1) =			% ktrans
%x(2) = 		% ve
x(3) = 1;		% unused
x(4) = resnorm;	% residual