%% Check the fit, as defined by multROI

function check = CHECKvpROI()
% Load the file

load('/home/data/scripts/multROI/ROIindices.mat');

% Choose the xdata file

[IMG,PathName,FilterIndex] = uigetfile([PathName '/*wpFIT*.mat'],'Choose fitted file');

load(fullfile(PathName, IMG));

load(fullfile(PathName, 'newAIF_with_vpFIT_voxels.mat'));
%size(x)
%roiindices
%save('temp.mat', 'roiindices', 'tumind');


ind = ismember(tumind, roiindices);
ind = find(ind == 1);

size(ind);

r = round(1 + (size(ind,1)-1).*rand(min(16,size(ind,1)),1))
a = figure;
clear y
for j = 1:numel(r)
    i = ind(r(j));
    timer = xdata{1}.timer;
    Ct = xdata{1}.Ct;
    Ct = Ct(:,i);
    
  
    %%FIX below
    y(j,:) = FXLStep1AIFcfit(x(i,1), x(i,2), xdata{1}.Cp, xdata{1}.timer)';
    
    
    figure(a)
    subplot(4,4,j), plot(timer, Ct, 'x'), hold on, plot(timer, y(j,:), 'r'), title(['Vo:' num2str(i), 'Ktve:' num2str(x(i,1:3)) 'R2:' num2str(x(i,4))])
    hold off
end

check = 1;