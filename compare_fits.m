
function compare_fits(results_d_path,background_image_path,show_original,show_ci)
load(results_d_path);
if xdata{1}.number_rois~=0
	compare_gui(xdata{1},fit_data,show_original,show_ci);
end
if fit_data.fit_voxels
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
% 	if ~exist('fitting_results','var')
% 		fit_data{1}.fitting_results = zeros(size(x,1),8);
% 		fit_data{1}.fitting_results(:,1:4) = x;
% 	end
% 	if ~exist('dce_model','var')
% 		dce_model = '';
% 	end

	while 1
		figure(1);
		[image_x,image_y]=ginput(1);
		image_x = int32(image_x);
		image_y = int32(image_y);

		size_y = size(background_image,1);
		% X and Y are reversed because of the transpose above
		image_id = image_x+(image_y-1)*size_y;
		voi = find(fit_data.tumind==image_id);
		if ~isempty(voi)
			plot_data.Ct = xdata{1}.Ct(:,voi);
			plot_data.Ct_original = xdata{1}.Ct_original(:,voi);
			plot_data.Cp = xdata{1}.Cp;
			plot_data.timer =  xdata{1}.timer;
			plot_data.fit_parameters = fit_data.fitting_results(voi,:);
			plot_data.dce_model = fit_data.dce_model;
			plot_data.show_original = show_original;
			plot_data.show_ci = show_ci;
			
			figure(2);
			plot_dce_curve(plot_data);
		end
	end
end

