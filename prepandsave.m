timeres = timereso;
xdata{1}.timer = timer;
xdata{1}.Ct    = Ct(starter:timeend,:);
numvoxels      = size(Ct,2);
%res = [0 0.25 0.25 2];
save(fullfile(PathName1, [rootname '_fitted_R1info.mat']));
fullfile(PathName1, [rootname '_fitted_R1info.mat'])