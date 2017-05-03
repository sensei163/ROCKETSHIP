function [x, residuals] = model_FXL_reference_region(Ct, Cp, timer, prefs)

options = fitoptions('Method', 'NonlinearLeastSquares',...
    'Algorithm', 'Trust-Region',...
    'MaxIter', prefs.MaxIter,...
    'MaxFunEvals', prefs.MaxFunEvals,...
    'TolFun', prefs.TolFun,...
    'TolX', prefs.TolX,...
    'Display', 'off',...
    'Lower',[prefs.lower_limit_ktrans prefs.lower_limit_ve prefs.lower_limit_ktrans_RR],...
    'Upper', [prefs.upper_limit_ktrans prefs.upper_limit_ve prefs.upper_limit_ktrans_RR],...
    'StartPoint', [prefs.initial_value_ktrans prefs.initial_value_ve prefs.initial_value_ktrans_RR],...
    'Robust', prefs.Robust);

ft = fittype('model_FXL_reference_region_cfit(Ktrans_TOI, ve_TOI, Ktrans_RR, Cp, T1)',...
    'independent', {'T1', 'Cp'},...
    'coefficients',{'Ktrans_TOI', 've_TOI','Ktrans_RR'});

[f, gof, output] = fit([timer, Cp'],Ct,ft, options);
confidence_interval = confint(f,0.95);

if output.exitflag<=0
    % Change start point to try for better fit
    new_options = fitoptions(options,...
        'StartPoint', [prefs.initial_value_ktrans*10 prefs.initial_value_ve prefs.initial_value_ktrans_RR]);
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
            'StartPoint', [prefs.initial_value_ktrans*100 prefs.initial_value_ve prefs.initial_value_ktrans_RR]);
        [new_f, new_gof, new_output] = fit([timer, Cp'],Ct,ft, new_options);

        if new_gof.sse < gof.sse
            f = new_f;
            gof = new_gof;
            output = new_output;
            confidence_interval = confint(f,0.95);
        end
    end
end

%Calculate the R2 fit
x(1) = f.Ktrans_TOI;		% ktrans_TOI
x(2) = f.ve_TOI;			% ve
x(3) = f.Ktrans_RR;         % ktrans_RR
x(4) = gof.sse;				% residual
x(5) = confidence_interval(1,1);% (95 lower CI of ktrans_TOI)
x(6) = confidence_interval(2,1);% (95 upper CI of ktrans_TOI)
x(7) = confidence_interval(1,2);% (95 lower CI of ve_TOI)
x(8) = confidence_interval(2,2);% (95 upper CI of ve_TOI)
x(9) = confidence_interval(1,3);% (95 lower CI of ktrans_RR)
x(10) = confidence_interval(2,3);% (95 upper CI of ktrans_RR)

residuals = output.residuals;

delete('tempve_RR.mat');