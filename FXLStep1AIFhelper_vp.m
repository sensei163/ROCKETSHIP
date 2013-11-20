%% Helper Function for nonlinear curvefit to FXLAIF model, vp
function x = FXLStep1AIFhelper_vp(xdata, voxel)
if nargin < 3
    smooth_model = 0;
end

% Fit the curve to brixfit
fitted = 0;

warning off
Ct = xdata{1}.Ct;
i  = voxel;
Ct = Ct(:,i);
if(fitted)
    [Ct bad] = smoothcurve(Ct, xdata, 0);
else
    bad = 1;
end

if smooth_model==1
	Ct = smooth(Ct,40,'moving');
elseif smooth_model==2
	Ct = smooth(Ct,0.05,'rlowess');
else
	% no smoothing
end

if(bad ~= -1)
	%configure the optimset for use with lsqcurvefit
	options = optimset('lsqcurvefit');

	%increase the number of function evaluations for more accuracy
	options.MaxFunEvals = 50;
	options.MaxIter     = 50;
	% options.TolFun      = 0.001;
	% options.TolX        = 0.001;
	options.TolFun      = 1e-8; %could use -8
	options.TolX        = 1e-4; %could use -4
	options.Diagnostics = 'off';
	options.Display     = 'off';
	options.Algorithm   = 'levenberg-marquardt';

	lb = [0 0 0];
	ub = [5 1 1];
% 	tic
	[x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@FXLStep1AIF_vp, ...
		[0.005 0.05 0.05], xdata, ...
		Ct',lb,ub,options);
% 	r2 = 1 - resnorm / norm(Ct-mean(Ct))^2;
	r2 = resnorm;
	
%	x(1) =			% ktrans
%	x(2) = 			% ve
% 	x(3) = 			% vp
	x(4) = resnorm;	% residual

% 	toc
else
    x(4) = -1;
end


