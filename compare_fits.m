function compare_fits(results_d_path,background_image_path,show_original,show_ci)

% compare_fits - Analyse the DCE fitting results from D_fit_voxels_func() 
%
% Inputs:
%  results_d_path        - *.mat Results from part D
%  background_image_path - path to image used to select voxels to see fit
%                          of
%  show_original         - show original unsmooted time curve
%  show_ci               - show approximate 95% confidence interval curve
%
% The script loads the data arrays generated from D_fit_voxels_func(). 
% Then one can visualize the fitted curves compared to the original time
% series data
% 
% Requires:
% compare_fits.m
% compare_gui.m
% plot_dce_curve.m
% FXLStep1AIF_vpcfit.m
% FXLStep1AIFcfit.m
% fxr_cfit.m
% niftitools
%
% 
% Samuel Barnes
% Caltech
% December 2013



load(results_d_path);
if fit_data.number_rois~=0
    compare_gui(xdata{1},fit_data,show_original,show_ci);
end
if fit_data.fit_voxels
    background_image = load_nii(background_image_path);
    background_image = double(background_image.img);

    if size(background_image,3)~=1
        background_image = background_image(:,:,floor(size(background_image,3)/2));
    end

    figure(1);
    % Note image transpose here
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
        % X and Y are reversed because of the transpose above
        [image_y,image_x]=ginput(1);
        image_y = int32(image_y);
        image_x = int32(image_x);

        size_y = size(background_image,1);
        image_id = image_y+(image_x-1)*size_y;
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
            plot_data.title = ['Voxel Location (' num2str(image_x) ',' num2str(image_y) ')'];
            
            if strcmp(fit_data.dce_model,'fxr')
                plot_data.R1o = xdata{1}.R1o(voi);
                plot_data.R1i = xdata{1}.R1i(voi);
                plot_data.r1 = xdata{1}.relaxivity;
                plot_data.fw = 0.8;
            end
            
            figure(2);
            plot_dce_curve(plot_data);
        end
    end
end

