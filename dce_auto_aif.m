% load('C:\Users\sbarnes\Documents\data\6 DCE Stroke\aging\Sam Analysis\Raw AIF\1 young - Copy\A-preAIFcalc.mat')

% TODO fix hard coded values
% TODO fix for 3D DYNAMIC

% Hard coded
end_ss = 110; %4.76 min
end_inject = 127; %5.5 min
t1_blood = 1.000; %secs
time_resolution = 2.6/60;
timer   = 0:time_resolution:time_resolution*(size(DYNAMIC,3)-1);
start_time = 1;
end_time = numel(timer);

% Shape data for detection
spatial_points = size(DYNAMIC,1)*size(DYNAMIC,2);
time_points = size(DYNAMIC,3);
dynam_pre_smooth = double(reshape(DYNAMIC,spatial_points,time_points)');
dynam_smooth = smooth(dynam_pre_smooth,25,'moving');
dynam_smooth = reshape(dynam_smooth,time_points,spatial_points);
% Scale from 0 to 1 so edge detector threshold is consistant
max_smooth = max(max(dynam_smooth));
min_smooth = min(min(dynam_smooth));
dynam_smooth = (dynam_smooth-min_smooth)./(max_smooth-min_smooth);

edge_detected = zeros(spatial_points,1);
all_edges = 0;
positive_edges = 0;
positive_edges_robust = 0;
good_fit_edges = 0;
tic
p = ProgressBar(100);
% Hack to make xdata compatible with parfor
xdata = cell(spatial_points,1);
parfor voxel_index = 1:spatial_points
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
            
            % Remove edges that fail after robust smoothing
            if max(is_edge_positive(end_ss:end_inject))>0
                positive_edges_robust = positive_edges_robust+1;

                % Fit voxel to biexponential and see what r^2 is
                aif_x = dynam_pre_smooth(:,voxel_index);
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
                if rsquare>0.8
                    edge_detected(voxel_index) = 1; 
                    good_fit_edges = good_fit_edges+1;                  
%                     figure;
%                     plot(timer,aif_x,'r.');
%                     hold on;
%                     plot(timer, aif_fitted,'b');
%                     title(num2str(rsquare));
%                     hold off;
%                     waitforbuttonpress
                end
            end
        end
    end
    if mod(voxel_index,100)==0, p.progress; end
end
p.stop;
clear xdata
toc
all_edges %#ok<NOPTS>
positive_edges %#ok<NOPTS>
positive_edges_robust %#ok<NOPTS>
good_fit_edges %#ok<NOPTS>

% figure(1)
% imshow(DYNAMIC(:,:,125)',[]);
% figure(2)
% imshow(reshape(edge_detected,size(DYNAMIC,1),size(DYNAMIC,2))');

% Plot found AIF voxels on image
figure;
imagesc(DYNAMIC(:,:,end_inject)');
colormap('gray');
size_image = size(DYNAMIC(:,:,end_inject)');
red_mask = cat(3, ones(size_image), zeros(size_image), zeros(size_image));
aif_mask = reshape(edge_detected,size(DYNAMIC,1),size(DYNAMIC,2));
hold on;
h_mask = imagesc(red_mask);
hold off;
set(h_mask, 'AlphaData',double(aif_mask'));



% Fit average AIF
aif_index = find(edge_detected>0);
if ~isempty(aif_index)
    aif_found = mean(dynam_pre_smooth(:,aif_index),2);
    aif_found = aif_found - mean(aif_found(1:end_ss));
    aif_found = aif_found./max(aif_found);
    xdata{1}.Cp    = aif_found;
    xdata{1}.timer = timer(start_time:end_time)';
    xdata{1}.step = [end_ss*time_resolution end_inject*time_resolution];
    % Show fit
    [aif_fitted, xAIF, xdataAIF, rsquare] = AIFbiexpfithelp(xdata, 1);
    figure;
    plot(timer,aif_found,'r.');
    hold on;
    plot(timer, aif_fitted,'b');
    title(num2str(rsquare));
    hold off;
else
    disp('No AIF found');
end
                

%%
% lvind = aif_ind;
% % T1LV     = (LV(lvind)+0.5).*max(Cp1000)';
% T1LV     = (LV(lvind)+0.5);
% T1       = T1LV;
% Stotallv = DYNAMLV;
% 
% 
% % Sss is the steady state Signal before injection
% Sss      = mean(Stotallv(round(steady_state_time(1)):round(steady_state_time(2)),:));
% Sstar    = ((1-exp(-tr./T1))./(1-cosd(fa).*exp(-tr./T1)));
% Stlv     = Stotallv;%(inj(2):end,:);
% 
% for j = 1:size(Stlv,2)
%     A(:,j) = 1-cosd(fa).*Sstar(j).*Stlv(:,j)./Sss(j);
%     B(:,j)  = 1-Sstar(j).*Stlv(:,j)./Sss(j);
% end
% 
% AB = A./B;
% % return
% % AB should not be less than 0. We interpolate the timeseries to clean
% % this. Threshold is 0.5;
% % up.
% [AB T1LV lvind BADspacelABv GOODspacelABv] = cleanAB(AB, T1LV,lvind, 'AIF', min(numel(T1LV)), 0.5);
% 
% R1tLV = double((1/tr).*log(AB));
% Sss = Sss(GOODspacelABv);
% Stlv = Stlv(:,GOODspacelABv);
% % R1 should be real. We interpolate the timeseries to clean
% % this. Threshold is 0.5;
% % up.
% [R1tLV T1LV lvind BADspacelv GOODspacelv] = cleanR1t(R1tLV, T1LV,lvind, 'AIF', min(numel(T1LV)), 0.5);
% Sss = Sss(GOODspacelv);
% Stlv = Stlv(:,GOODspacelv);
% 
% for j = 1:numel(T1LV)
%     % Scale Sss values to T1 values
%     ScaleFactorlv = (1/T1LV(j)) - mean(R1tLV(round(steady_state_time(1)):round(steady_state_time(2)), j));
%     R1tLV(:,j) = R1tLV(:,j) + ScaleFactorlv;
% end
% 
% 
% % 9. Convert to concentrations
% 
% % f.1 AIF
% %if(~iterremove)
% Cp = (R1tLV-repmat((1./T1LV)', [size(R1tLV, 1) 1]))./(relaxivity*(1-hematocrit));