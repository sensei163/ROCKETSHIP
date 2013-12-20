%% Check the fit

function check = CHECKFXRROI()

load('/home/data/scripts/multROI/ROIindices.mat');

% Choose the xdata file

[IMG,PathName,FilterIndex] = uigetfile([PathName '/*FXR*.mat'],'Choose fitted file');

load(fullfile(PathName, IMG));

load(fullfile(PathName, 'newAIF_FXR_voxels.mat'));

ind = ismember(tumind, roiindices);
ind = find(ind == 1);

r = round(1 + (size(ind,1)-1).*rand(min(16,size(ind,1)),1));
a = figure;
clear y
for j = 1:numel(r)
    i = ind(r(j));
    timer = xdata{1}.timer;
    Rt1 = xdata{1}.R1tTOI;
    Rt1 = Rt1(:,i);
    
    y(j,:) =  FXRAIF(x(i,:), xdata); %FXLStep1AIFcfit(x(i,1), x(i,2), xdata{1}.Cp, xdata{1}.timer);
    
    
    figure(a)
    subplot(4,4,j), plot(timer, Rt1, 'x'), hold on, plot(timer, y(j,:), 'r'), title(['Voxel: ' num2str(i), 'K_{trans} v_{e} /tau_{i}: ' num2str(x(i,1:3)) 'R^2 fit: ' num2str(x(i,end))])
    hold off
end

check = 1;