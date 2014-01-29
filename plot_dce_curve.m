function plot_dce_curve(plot_data)
fit_parameters = plot_data.fit_parameters;

% fit_parameters = fitting_results(voi,:);
if strcmp(plot_data.model_name,'aif')
    fit_curve = FXLStep1AIFcfit(fit_parameters(1),fit_parameters(2),plot_data.Cp,plot_data.timer);
    fit_curve_low = FXLStep1AIFcfit(fit_parameters(5),fit_parameters(7),plot_data.Cp,plot_data.timer);
    fit_curve_high = FXLStep1AIFcfit(fit_parameters(6),fit_parameters(8),plot_data.Cp,plot_data.timer);
elseif strcmp(plot_data.model_name,'aif_vp')
    fit_curve = FXLStep1AIF_vpcfit(fit_parameters(1),fit_parameters(2),fit_parameters(3),plot_data.Cp,plot_data.timer);
    fit_curve_low = FXLStep1AIF_vpcfit(fit_parameters(5),fit_parameters(7),fit_parameters(9),plot_data.Cp,plot_data.timer);
    fit_curve_high = FXLStep1AIF_vpcfit(fit_parameters(6),fit_parameters(8),fit_parameters(10),plot_data.Cp,plot_data.timer);
elseif strcmp(plot_data.model_name,'fxr')
    fit_curve = fxr_cfit(fit_parameters(1),fit_parameters(2),fit_parameters(3),plot_data.Cp,plot_data.timer,plot_data.R1o,plot_data.R1i,plot_data.r1,plot_data.fw);
    fit_curve_low = fxr_cfit(fit_parameters(5),fit_parameters(7),fit_parameters(9),plot_data.Cp,plot_data.timer,plot_data.R1o,plot_data.R1i,plot_data.r1,plot_data.fw);
    fit_curve_high = fxr_cfit(fit_parameters(6),fit_parameters(8),fit_parameters(10),plot_data.Cp,plot_data.timer,plot_data.R1o,plot_data.R1i,plot_data.r1,plot_data.fw);
else
    fit_curve = FXLStep1AIFcfit(fit_parameters(1),fit_parameters(2),plot_data.Cp,plot_data.timer);
    fit_curve_low = FXLStep1AIFcfit(fit_parameters(5),fit_parameters(7),plot_data.Cp,plot_data.timer);
    fit_curve_high = FXLStep1AIFcfit(fit_parameters(6),fit_parameters(8),plot_data.Cp,plot_data.timer);
end

if plot_data.show_original
    p_original = plot(plot_data.timer,plot_data.Ct_original);
    set(p_original(1),'Color',[0.7 0.7 0.8]);
    hold on
end

p = plot(plot_data.timer,plot_data.Ct,...
    plot_data.timer,fit_curve(:,1));
% 		set(p(1),'Marker','.','MarkerSize',4,'LineStyle','none');
set(p(1),'Color',[0.6 0.6 0.7]);
set(p(2),'Color','b');
if plot_data.show_original
    set(p(1),'Color','k');
end
if plot_data.show_ci
    hold on
    p_ci = plot(plot_data.timer,fit_curve_low(:,1),...
        plot_data.timer,fit_curve_high(:,1));
    set(p_ci(1),'Color','k');
    set(p_ci(2),'Color','k');
end
hold off


plot_limits = axis;
% Estimate Error
if strcmp(plot_data.model_name,'aif')
    ktrans_error = mean([abs(fit_parameters(1)-fit_parameters(5)) abs(fit_parameters(1)-fit_parameters(6))]);
    ve_error = mean([abs(fit_parameters(2)-fit_parameters(7)) abs(fit_parameters(2)-fit_parameters(8))]);
    plot_str(1) = {[' Ktrans = ' num2str(fit_parameters(1),2) '±' num2str(ktrans_error,2)]};
    plot_str(2) = {[' Ve = ' num2str(fit_parameters(2),2) '±' num2str(ve_error,2)]};
    plot_str(3) = {[' residual = ' num2str(fit_parameters(4))]};
elseif strcmp(plot_data.model_name,'aif_vp')
    ktrans_error = mean([abs(fit_parameters(1)-fit_parameters(5)) abs(fit_parameters(1)-fit_parameters(6))]);
    ve_error = mean([abs(fit_parameters(2)-fit_parameters(7)) abs(fit_parameters(2)-fit_parameters(8))]);
    vp_error = mean([abs(fit_parameters(3)-fit_parameters(9)) abs(fit_parameters(3)-fit_parameters(10))]);
    plot_str(1) = {[' Ktrans = ' num2str(fit_parameters(1),2) '±' num2str(ktrans_error,2)]};
    plot_str(2) = {[' Ve = ' num2str(fit_parameters(2),2) '±' num2str(ve_error,2)]};
    plot_str(3) = {[' Vp = ' num2str(fit_parameters(3),2) '±' num2str(vp_error,2)]};
    plot_str(4) = {[' residual = ' num2str(fit_parameters(4))]};
elseif strcmp(plot_data.model_name,'fxr')
    ktrans_error = mean([abs(fit_parameters(1)-fit_parameters(5)) abs(fit_parameters(1)-fit_parameters(6))]);
    ve_error = mean([abs(fit_parameters(2)-fit_parameters(7)) abs(fit_parameters(2)-fit_parameters(8))]);
    vp_error = mean([abs(fit_parameters(3)-fit_parameters(9)) abs(fit_parameters(3)-fit_parameters(10))]);
    plot_str(1) = {[' Ktrans = ' num2str(fit_parameters(1),2) '±' num2str(ktrans_error,2)]};
    plot_str(2) = {[' Ve = ' num2str(fit_parameters(2),2) '±' num2str(ve_error,2)]};
    plot_str(3) = {[' tau = ' num2str(fit_parameters(3),2) '±' num2str(vp_error,2)]};
    plot_str(4) = {[' residual = ' num2str(fit_parameters(4))]};
else
    ktrans_error = mean([abs(fit_parameters(1)-fit_parameters(5)) abs(fit_parameters(1)-fit_parameters(6))]);
    ve_error = mean([abs(fit_parameters(2)-fit_parameters(7)) abs(fit_parameters(2)-fit_parameters(8))]);
    plot_str(1) = {[' Ktrans = ' num2str(fit_parameters(1),2) '±' num2str(ktrans_error,2)]};
    plot_str(2) = {[' Ve = ' num2str(fit_parameters(2),2) '±' num2str(ve_error,2)]};
    plot_str(3) = {[' residual = ' num2str(fit_parameters(4))]};
end
text(plot_limits(1),plot_limits(3),plot_str,...
    'Color', 'black',...
    'VerticalAlignment','bottom',...
    'HorizontalAlignment','left');
title(plot_data.title);
xlabel(plot_data.x_units);
ylabel(plot_data.y_units);