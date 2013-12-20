%% 3. Fit to model VP on a voxel by voxel basis

function xdata = D_fit_SAUC_voxels_COH(directory1, CPU, CHECK, r2filter, timestamp)
%% 1. Load files from C_fit_AIFnoCP.m

% r2filter = 0; % Filter out all fits with r2 < r2filter
% CPU      = 0; % Use only if you are doing multicore
CPU

% Load the data files

rootname  = 'SAUC';
% [gogo,PathName1,FilterIndex] = uigetfile(['/data/studies/' '/*AIF_with_vpFIT_ROI.mat'],'Choose R1 file');
%
% % Load the data files
% directory = PathName1
% rootname  = strrep(gogo, '.nii', '');
% load(fullfile(PathName1, gogo));
directory1 = strrep(directory1, 'AIF_with_vpFIT_ROI.mat', '_fitted_R1info.mat');
load(directory1)

inject = 2.5;
timer = [0:2/60:25];

timer = timer(1:size(Stotallv,1));
ind1    = find(timer >= inject);
ind1  = ind1(1);
ind2    = find(timer<=timestamp);
ind2  = ind2(end);
% AUC AIF
AIF = mean(Stotallv,2);
AIF = AIF-mean(Sss);

AIFAUC = trapz(timer(ind1:ind2), AIF(ind1:ind2));

Stotaltum = Stotaltum(:, 1:numel(tumind));


for i = 1:size(Stotaltum,2)
    xdata{1}.timer = timer;
    xdata{1}.inject = inject;
    
    
    CurCt = Stotaltum(:,i);
    
    % [CurCt bad] = smoothcurve(CurCt, xdata, 1);
    CurCt = CurCt-mean(Ssstum(i));
    
    DiffY = diff(CurCt);
    DiffX = diff(timer);
    
    Diff(i) = max(DiffY(:)./DiffX(:));
    
    
    bad = 1;
    if(bad ~= -1)
        AUC(i) = trapz(timer(ind1:ind2), CurCt(ind1:ind2));
    else
        AUC(i) = -1;
        Diff(i) = -1;
    end
    
end


badind = find(AUC<0);
AUC(badind) = [];
tumind(badind) = [];
Diff(badind) = [];


%Normalize AUC
NAUC = AUC./AIFAUC;

AIFAUC
%% f) Make maps

AUCROI    = zeros(size(currentimg));

DiffROI    = zeros(size(currentimg));



AUCROI(tumind) = NAUC;
DiffROI(tumind)= Diff;

xdata{1}.AUC = AUC;
xdata{1}.NAUC= NAUC;
xdata{1}.Diff= Diff;
xdata{1}.AUCCp=AIF;
xdata{1}.tumind = tumind;


res = [0 0.25 0.25 2];

%% Calculate Fractals
% Take central slice

A = AUCROI(:,:,2);
[i j] = find(A~=0);
BLOCKA = zeros((max(i)-min(i)+1), max(j)-min(j)+1);
BLOCKA(1:end, 1:end) = A(min(i):max(i), min(j):max(j));

%imagesc(A)
fract = FractalDim(BLOCKA,10);


xdata{1}.fract = fract;

%% h) Save image files

[discard actual] = fileparts(strrep(dynamname.fileprefix, '\', '/'));

save_nii(make_nii(AUCROI, res(2:4), [1 1 1]), fullfile(fileparts(directory1), [actual '_SAUC_min' num2str(timestamp) '.nii']));
save_nii(make_nii(DiffROI, res(2:4), [1 1 1]), fullfile(fileparts(directory1), [actual '_SDiff_min' num2str(timestamp) '.nii']));

a = 1;





