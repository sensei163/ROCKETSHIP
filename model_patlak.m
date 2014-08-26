function [x, residuals] = model_patlak(Ct,Cp,timer,prefs)

% Use Curvefitting tool box instead of optimization toolbox (lsqcurvefit)
% as curvefitting will easily return confidence intervals on the fit
% performance of the two appears to be the same
options = fitoptions('Method', 'NonlinearLeastSquares',...
    'Algorithm', 'Trust-Region',...
    'MaxIter', prefs.MaxIter,...
    'MaxFunEvals', prefs.MaxFunEvals,...
    'TolFun', prefs.TolFun,...
    'TolX', prefs.TolX,...
    'Display', 'off',...
    'Lower',[prefs.lower_limit_ktrans prefs.lower_limit_vp],...
    'Upper', [prefs.upper_limit_ktrans prefs.upper_limit_vp],...
    'StartPoint', [prefs.initial_value_ktrans prefs.initial_value_vp],...
    'Robust', prefs.Robust);
ft = fittype('model_patlak_cfit( Ktrans, vp, Cp, time)',...
    'independent', {'time', 'Cp'},...
    'coefficients',{'Ktrans', 'vp'});
[f, gof, output] = fit([timer, Cp'],Ct,ft, options);
confidence_interval = confint(f,0.95);

if output.exitflag<=0
    % Change start point to try for better fit
    new_options = fitoptions(options,...
        'StartPoint', [prefs.initial_value_ktrans*10 prefs.initial_value_vp]);
    [new_f, new_gof, new_output] = fit([timer, Cp'],Ct,ft, new_options);
    
    if new_gof.sse < gof.sse
        f = new_f;
        gof = new_gof;
        output = new_output;
        confidence_interval = confint(f,0.95);
    end
    
    if output.exitflag<=0
        % Change start point to try for better fit
        new_options = fitoptions(options,...
            'StartPoint', [prefs.initial_value_ktrans*100 prefs.initial_value_vp]);
        [new_f, new_gof, new_output] = fit([timer, Cp'],Ct,ft, new_options);

        if new_gof.sse < gof.sse
            f = new_f;
            gof = new_gof;
            output = new_output;
            confidence_interval = confint(f,0.95);
        end
    end
end

x(1) = f.Ktrans;			% ktrans
x(2) = f.vp;				% vp
x(3) = gof.sse;				% residual
x(4) = confidence_interval(1,1);% (95 lower CI of ktrans)
x(5) = confidence_interval(2,1);% (95 upper CI of ktrans)
x(6) = confidence_interval(1,2);% (95 lower CI of vp)
x(7) = confidence_interval(2,2);% (95 upper CI of vp)

residuals = output.residuals;
