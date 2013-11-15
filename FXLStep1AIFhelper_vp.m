%% Helper Function for nonlinear curvefit to FXLAIF model, vp
function x = FXLStep1AIFhelper_vp(xdata, voxel)

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

if(bad ~= -1)
	%configure the optimset for use with lsqcurvefit
	options = optimset('lsqcurvefit');

	%increase the number of function evaluations for more accuracy
	options.MaxFunEvals = 50;
	options.MaxIter     = 50;
	% options.TolFun      = 0.001;
	% options.TolX        = 0.001;
	options.TolFun      = 1e-10; %could use -8
	options.TolX        = 1e-5; %could use -4
	options.Diagnostics = 'off';
	options.Display     = 'off';
	options.Algorithm   = 'levenberg-marquardt';

	lb = [0 0 0];
	ub = [5 1 1];
	%
% 	tic
	%try
	[x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@FXLStep1AIF_vp, ...
		[0.005 0.05 0.05], xdata, ...
		Ct','','',options);
% 	r2 = 1 - resnorm / norm(Ct-mean(Ct))^2;
	r2 = resnorm;
	x(end+1) = r2;

	% matlabpool local 7
	% RS = RandomStartPointSet('NumStartPoints', 50);
	% options = optimset('Algorithm','trust-region-reflective');
	% problem = createOptimProblem('lsqcurvefit', 'x0', x(1:3), 'objective',@FXLStep1AIF_vp, 'lb', [0 0 0], 'ub', [1 1 1], 'xdata', xdata, 'options', options, 'ydata', Ct' );
	% ms = MultiStart('StartPointsToRun', 'bounds', 'UseParallel', 'always', 'Display', 'off');
	% [xming,fming,flagg,outptg,manyminsg] = run(ms,problem, RS);
	% disp('Applying Multistart Global optimization KtransRR, Ktrans TOI, ve TOI:'),
	% x(1:3) = xming
	% matlabpool close

	% Use Curvefitting tool box

	%options = fitoptions('Method', 'NonlinearLeastSquares', 'Algorithm', 'Levenberg-Marquardt', 'MaxIter', 150, 'MaxFunEvals', 150, 'TolFun', 0.001, 'TolX', 0.001, 'Display', 'off', 'Lower',[0 0 0], 'Upper', [5 1 1], 'StartPoint', [0.5 0.5 0.5]);


	%ft = fittype('FXLStep1AIF_vpcfit( Ktrans, ve, vp, Cp, T1)', 'independent', {'T1', 'Cp'}, 'coefficients',{'Ktrans', 've', 'vp'});
	%{
		[f gof] = fit([xdata{1}.timer', xdata{1}.Cp],Ct,ft, options);


		%Calculate the R2 fit
	x(1) = f.Ktrans;
	x(2) = f.ve;
	x(3) = f.vp;
	x(end+1) = gof.rsquare;
	x(end+1) = gof.adjrsquare;

	%}
	%catch
	%x = [-1 -1 i i];
	%end
% 	toc
else
    x(4) = -1;
end


