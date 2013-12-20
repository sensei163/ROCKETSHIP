%% Helper Function for nonlinear curvefit to FXR model, no vp
function x = FXRAIFhelper(xdata, voxel)
warning off
    fitted = 0;
R1t = xdata{1}.R1t;
i  = voxel;
R1t = R1t(:,i);

%figure, plot(R1t)
if(fitted)
[R1t bad] = smoothcurve(R1t, xdata, 1);
else
    bad = 1;
end

if(bad ~= -1)
%configure the optimset for use with lsqcurvefit
options = optimset('lsqcurvefit');

%increase the number of function evaluations for more accuracy
options.MaxFunEvals = 50;
options.MaxIter     = 50;
options.TolFun      = 0.001;
options.TolX        = 0.001;
options.Diagnostics = 'off';
options.Display     = 'off';
options.Algorithm   = 'levenberg-marquardt';

lb = [0 0 0];
ub = [5 1 5];

%tic 
try
[x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@FXRAIF, ...
    [0.5 0.5 0.5], xdata, ...
    R1t,'','',options);
r2 = 1 - resnorm / norm(R1t-mean(R1t))^2;
x(end+1) = r2;

catch
    x = [-1 -1 i i];
end
% Use Curvefitting tool box
% 
% options = fitoptions('Method', 'NonlinearLeastSquares', 'Algorithm', 'Levenberg-Marquardt', 'MaxIter', 150, 'MaxFunEvals', 150, 'TolFun', 0.001, 'TolX', 0.001, 'Display', 'off', 'Lower',[0 0 0], 'Upper', [5 1 5], 'StartPoint', [0.5 0.5 0.5]);
% 
% 
% ft = fittype('FXRAIF(Ktrans, ve, tau1, R1o, R1i, fw, r1, Cp, T1)', 'independent', {'T1', 'Cp'},'problem', {'R1o', 'R1i', 'fw', 'r1'}, 'coefficients',{'Ktrans', 've', 'tau1'});
% [f gof] = fit([ xdata{1}.timer', xdata{1}.Cp],R1t,ft,'problem', {R1o, R1i, fw, r1}, options)
%toc

%Calculate the R2 fit
% x(1) = f.Ktrans;
% x(2) = f.ve;
% x(3) = f.tau1;
% x(end+1) = gof.rsquare;
% x(end+1) = gof.adjrsquare;

else
    x(4) = -1;
end
