
function compare_fits(results_d_path,background_image_path,show_original,show_ci)

load(results_d_path);
background_image = load_nii(background_image_path);
background_image = double(background_image.img);

if size(background_image,3)~=1
	background_image = background_image(:,:,floor(size(background_image,3)/2));
end

figure(1);
imshow(background_image',...
	[prctile(reshape(background_image,1,[]),5) prctile(reshape(background_image,1,[]),95)],...
	'InitialMagnification', 400);

% For backwards compatability with old results
if ~exist('fitting_results','var')
	fitting_results = zeros(size(x,1),8);
	fitting_results(:,1:4) = x;
end
if ~exist('dce_model','var')
	dce_model = '';
end

while 1
	figure(1);
	[image_x,image_y]=ginput(1);
	image_x = int32(image_x);
	image_y = int32(image_y);
	
	size_y = size(background_image,1);
	% X and Y are reversed because of the transpose above
	image_id = image_x+(image_y-1)*size_y;
	voi = find(tumind==image_id);
	if ~isempty(voi)
		Ct = xdata{1,1}.Ct(:,voi);

% 		if exist('time_smoothing','var')
% 			if strcmp(time_smoothing,'moving')
% 				Ct = smooth(Ct,time_smoothing_window,'moving');
% 			elseif strcmp(time_smoothing,'rlowess')
% 				Ct = smooth(Ct,time_smoothing_window/size(Ct,1),'rlowess');
% 			else
% 				% no smoothing
% 			end
% 		else
% 			% Defaut smoothing
% 			Ct = smooth(Ct,9/size(Ct,1),'rlowess');
% 		end

		figure(2);
% 		fit_parameters = FXLStep1AIFhelper(xdata,voi,0);
		fit_parameters = fitting_results(voi,:);
		if strcmp(dce_model,'aif')
			fit_curve = FXLStep1AIFcfit(fit_parameters(1),fit_parameters(2),xdata{1}.Cp,xdata{1}.timer);
			fit_curve_low = FXLStep1AIFcfit(fit_parameters(5),fit_parameters(7),xdata{1}.Cp,xdata{1}.timer);
			fit_curve_high = FXLStep1AIFcfit(fit_parameters(6),fit_parameters(8),xdata{1}.Cp,xdata{1}.timer);
		elseif strcmp(dce_model,'aif_vp')
			fit_curve = FXLStep1AIF_vpcfit(fit_parameters(1),fit_parameters(2),fit_parameters(3),xdata{1}.Cp,xdata{1}.timer);
			fit_curve_low = FXLStep1AIF_vpcfit(fit_parameters(5),fit_parameters(7),fit_parameters(9),xdata{1}.Cp,xdata{1}.timer);
			fit_curve_high = FXLStep1AIF_vpcfit(fit_parameters(6),fit_parameters(8),fit_parameters(10),xdata{1}.Cp,xdata{1}.timer);
		else
			fit_curve = FXLStep1AIFcfit(fit_parameters(1),fit_parameters(2),xdata{1}.Cp,xdata{1}.timer);
			fit_curve_low = FXLStep1AIFcfit(fit_parameters(5),fit_parameters(7),xdata{1}.Cp,xdata{1}.timer);
			fit_curve_high = FXLStep1AIFcfit(fit_parameters(6),fit_parameters(8),xdata{1}.Cp,xdata{1}.timer);
		end
		
		if show_original
			p_original = plot(xdata{1}.timer,xdata{1}.Ct_original(:,voi));
			set(p_original(1),'Color',[0.7 0.7 0.8]);
			hold on
		end
		
		p = plot(xdata{1}.timer,Ct,...
			xdata{1}.timer,fit_curve(:,1));
% 		set(p(1),'Marker','.','MarkerSize',4,'LineStyle','none');
		set(p(1),'Color',[0.6 0.6 0.7]);
		set(p(2),'Color','b');
		if show_original
			set(p(1),'Color','k');
		end
		if show_ci
			hold on
			p_ci = plot(xdata{1}.timer,fit_curve_low(:,1),...
				xdata{1}.timer,fit_curve_high(:,1));
			set(p_ci(1),'Color','k');
			set(p_ci(2),'Color','k');
		end
		hold off
		

		plot_limits = axis;
		% Estimate Error
		if strcmp(dce_model,'aif')
			ktrans_error = mean([abs(fit_parameters(1)-fit_parameters(5)) abs(fit_parameters(1)-fit_parameters(6))]);
			ve_error = mean([abs(fit_parameters(2)-fit_parameters(7)) abs(fit_parameters(2)-fit_parameters(8))]);
			plot_str(1) = {[' Ktrans = ' num2str(fit_parameters(1),2) '±' num2str(ktrans_error,2)]};
			plot_str(2) = {[' Ve = ' num2str(fit_parameters(2),2) '±' num2str(ve_error,2)]};
			plot_str(3) = {[' residual = ' num2str(fit_parameters(4))]};
		elseif strcmp(dce_model,'aif_vp')
			ktrans_error = mean([abs(fit_parameters(1)-fit_parameters(5)) abs(fit_parameters(1)-fit_parameters(6))]);
			ve_error = mean([abs(fit_parameters(2)-fit_parameters(7)) abs(fit_parameters(2)-fit_parameters(8))]);
			vp_error = mean([abs(fit_parameters(3)-fit_parameters(9)) abs(fit_parameters(3)-fit_parameters(10))]);
			plot_str(1) = {[' Ktrans = ' num2str(fit_parameters(1),2) '±' num2str(ktrans_error,2)]};
			plot_str(2) = {[' Ve = ' num2str(fit_parameters(2),2) '±' num2str(ve_error,2)]};
			plot_str(3) = {[' Vp = ' num2str(fit_parameters(3),2) '±' num2str(vp_error,2)]};
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
	end
end

