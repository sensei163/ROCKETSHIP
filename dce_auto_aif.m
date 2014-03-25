function [aif_index, dynamic_output] = dce_auto_aif(dynamic_input, mask_index,dimx,dimy,dimz)
% load('C:\Users\sbarnes\Documents\data\6 DCE Stroke\aging\Sam Analysis\Raw AIF\1 young - Copy\A-preAIFcalc.mat')
% load('C:\Users\sbarnes\Documents\data\6 DCE Stroke\aging\Sam Analysis\Fitted Auto AIF\1 young\A-preAIFcalcPostDrift.mat')
% load('C:\Users\sbarnes\Documents\data\6 DCE Stroke\aging\Sam Analysis\Fitted Auto AIF\2 old\A-preAutoFitWithDrift.mat')

% TODO fix hard coded values
% TODO fix for 3D DYNAMIC
disp('Running auto AIF selection, searching for AIF...');
% Hard coded
end_ss = 110; %4.76 min
end_inject = 127; %5.5 min
time_resolution = 2.6/60;
timer   = 0:time_resolution:time_resolution*(size(dynamic_input,1)-1);
start_time = 1;
end_time = numel(timer);

% Shape data for detection
spatial_points = size(dynamic_input,2);
time_points = size(dynamic_input,1);
% dynam_pre_smooth = double(reshape(DYNAMIC,spatial_points,time_points)');
% dynamic_input = dynamic_input;
dynam_smooth = smooth(dynamic_input,25,'moving');
dynam_smooth = reshape(dynam_smooth,time_points,spatial_points);
% Scale from 0 to 1 so edge detector threshold is consistant
max_smooth = max(max(dynam_smooth));
min_smooth = min(min(dynam_smooth));
dynam_smooth = (dynam_smooth-min_smooth)./(max_smooth-min_smooth);

dynamic_output = [];
edge_detected = zeros(spatial_points,1);
all_edges = 0;
positive_edges = 0;
% positive_edges_robust = 0;
good_fit_edges = 0;
tic
p = ProgressBar(100);
% Hack to make xdata compatible with parfor
xdata = cell(spatial_points,1);
parfor voxel_index = 1:spatial_points
    % Only consider points that are inside the mask
%     if isempty(find(voxel_index==mask_index,1))
%         continue;
%     end
    
    % Edge detect with sobel filter
    [bw thresh]= edge(dynam_smooth(:,voxel_index),'sobel',0.007);
    % Look only for edges within injection window
    if max(bw(end_ss:end_inject))==1   
        all_edges = all_edges+1;
        
        
%         bw_plot = bw*0.5;     
% %         bw_plot = is_edge_positive_robust.*10+min(dynam_smooth(:,voxel_index));
%         figure(1); 
%         plot(dynam_smooth(:,voxel_index)); hold on
% %         plot(dynam_smooth_robust); hold on
%         plot(bw_plot,'r');
%         hold off;
% %         figure(2);
% %         plot(dynam_pre_smooth(:,voxel_index));
% %         thresh
% %         waitforbuttonpress
        
        
        % Calculate gradient to see if edge is positive
        op = fspecial('sobel')/8;
        x_mask = op;
        bx = imfilter(dynam_smooth(:,voxel_index),x_mask,'replicate');
        is_edge_positive = -bw.*bx;

        % Eliminate negative edges
        if max(is_edge_positive(end_ss:end_inject))>0
            positive_edges = positive_edges+1;
            % Run robust smoothing to remove more noise
            % Very slow so only do it after is has been pre selected
            % with faster algorithm
%             dynam_smooth_robust = smooth(dynam_pre_smooth(:,voxel_index),25,'rlowess');
%             max_smooth = max(dynam_smooth_robust);
%             min_smooth = min(dynam_smooth_robust);
%             dynam_smooth_robust = (dynam_smooth_robust-min_smooth)./(max_smooth-min_smooth);
% 
%             bw_robust = edge(dynam_smooth_robust,'sobel',0.007);
%             is_edge_positive_robust = -bw_robust.*bx;
%             
%             % Remove edges that fail after robust smoothing
%             if max(is_edge_positive_robust(end_ss:end_inject))>0
%                 positive_edges_robust = positive_edges_robust+1;

            % Fit voxel to biexponential and see what r^2 is
            aif_x = dynamic_input(:,voxel_index);
            % Set baseline to zero so it will fit properly
            aif_x = aif_x - mean(aif_x(1:end_ss));
            % Scale to 1
            aif_x = aif_x./max(aif_x);
            xdata{voxel_index}.Cp    = aif_x;
            xdata{voxel_index}.timer = timer(start_time:end_time)';
            xdata{voxel_index}.step = [end_ss*time_resolution end_inject*time_resolution];
            % Run fit
            [aif_fitted, ~, ~, rsquare] = AIFbiexpfithelp(xdata{voxel_index}, 0);
            % Remove bad fits
            if rsquare>0.8 && aif_fitted(end)<0.3
                edge_detected(voxel_index) = 1; 
                good_fit_edges = good_fit_edges+1;                  
%                 figure;
%                 plot(timer,aif_x,'r.');
%                 hold on;
%                 plot(timer, aif_fitted,'b');
%                 title(num2str(rsquare));
%                 hold off;
%                 waitforbuttonpress
            end
%             end
        end
    end
    if mod(voxel_index,100)==0, p.progress; end
end
p.stop;
% Free memeory
clear xdata
disp('Finished AIF search');
disp(['Found ' num2str(good_fit_edges) ' AIF voxels']);
toc
% all_edges %#ok<NOPRT>
% positive_edges %#ok<NOPRT>
% positive_edges_robust %#ok<NOPRT>
% good_fit_edges %#ok<NOPRT>

% figure(1)
% imshow(DYNAMIC(:,:,125)',[]);
% figure(2)
% imshow(reshape(edge_detected,size(DYNAMIC,1),size(DYNAMIC,2))');

% Plot found AIF voxels on image
figure;
dynamic_slice = zeros(dimx,dimy);
dynamic_slice(mask_index) = dynamic_input(end_inject,:);
imagesc(dynamic_slice'), axis off;
colormap('gray');
red_mask = cat(3, ones([dimy dimx]), zeros([dimy dimx]), zeros([dimy dimx]));
aif_mask = zeros(dimx,dimy);
aif_mask(mask_index) = edge_detected;
% aif_mask = reshape(edge_detected,dimx,dimy);
hold on;
h_mask = imagesc(red_mask);
hold off;
set(h_mask, 'AlphaData',double(aif_mask'));



% Fit average AIF
aif_index = find(aif_mask>0);
dynamic_output = dynamic_input(:,edge_detected>0);

if ~isempty(aif_index)
    aif_found = mean(dynamic_output,2);
    aif_found = aif_found - mean(aif_found(1:end_ss));
    aif_found = aif_found./max(aif_found);
    xdata{1}.Cp    = aif_found;
    xdata{1}.timer = timer(start_time:end_time)';
    xdata{1}.step = [end_ss*time_resolution end_inject*time_resolution];
    % Show fit
    [aif_fitted, xAIF, xdataAIF, rsquare] = AIFbiexpfithelp(xdata, 0);
    figure;
    plot(timer,aif_found,'r.');
    hold on;
    plot(timer, aif_fitted,'b');
    title(['r^2 = ' num2str(rsquare)]);
    hold off;
else
    disp('No AIF found');
end
       

