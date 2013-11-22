
function compare_fits(results_d_path,background_image_path)

% base_directory = 'C:\Users\sbarnes\Documents\data\6 DCE Stroke\sb01_06nov13.mH1';
% load('C:\Users\sbarnes\Documents\data\6 DCE Stroke\sb01_06nov13.mH1\20131106_000_sb01_06nov13_dynamic_01_s_AIF_with_vpFIT_ROI.mat')
% % load('C:\Users\sbarnes\Documents\data\6 DCE Stroke\sb01_06nov13.mH1\r20\20131106_000_sb01_06nov13_dynamic_01_s_AIF_with_vpFIT_ROI.mataif_FIT_voxels.mat')
% [gogo,PathName,FilterIndex] = uigetfile([base_directory '/*aif_FIT_voxels'],'Choose results file');
% if gogo==0
% 	disp('User selected cancel');
% 	return
% end
% load(fullfile(PathName, gogo));

load(results_d_path);

% place = '';
% [tumor,PathName3,FilterIndex] = uigetfile([base_directory '/*.nii'],'Choose Background Image');
% if tumor==0
% 	return
% end
% background_image = load_nii(fullfile(PathName3, place, tumor));
background_image = load_nii(background_image_path);
background_image = double(background_image.img);

% background_image = zeros(size(currentimg));
% background_image(tumind) = x(:,1);
if size(background_image,3)~=1
	background_image = background_image(:,:,floor(size(background_image,3)/2));
end

figure(1);
% imshow(background_image'.*1000,'InitialMagnification', 400);
imshow(background_image',...
	[prctile(reshape(background_image,1,[]),5) prctile(reshape(background_image,1,[]),95)],...
	'InitialMagnification', 400);

while 1
	figure(1);
	[image_x,image_y]=ginput(1);
	image_x = int32(image_x);
	image_y = int32(image_y);
	
	size_y = size(background_image,1);
	% X and Y are reversed because of the transpose above
% 	image_id = image_y+(image_x-1)*size_y;
	image_id = image_x+(image_y-1)*size_y;
	voi = find(tumind==image_id);
	if ~isempty(voi)
		curve = xdata{1,1}.Ct(:,voi);
% 		moving_curve = smooth(curve,8,'moving');
		lowess_curve = smooth(curve,0.01,'rlowess');

		figure(2);
% 		fit_parameters = FXLStep1AIFhelper(xdata,voi,0);
		fit_parameters = x(voi,:);
		fit_curve = FXLStep1AIF(fit_parameters(1:3),xdata);
		p = plot(1:size(lowess_curve,1),lowess_curve,1:size(fit_curve,2),fit_curve(1,:));
		plot_limits = axis;
		plot_str(1) = {[' Ktrans = ' num2str(fit_parameters(1))]};
		plot_str(2) = {[' Ve = ' num2str(fit_parameters(2))]};
		plot_str(3) = {[' residual = ' num2str(fit_parameters(4))]};
		text(plot_limits(1),plot_limits(3),plot_str,...
			'Color', 'black',...
			'VerticalAlignment','bottom',...
			'HorizontalAlignment','left');
	end
end

