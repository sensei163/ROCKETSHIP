function results = A_make_R1maps_func(dce_path,t1_aif_path,t1_roi_path,noise_path,tr,fa,hematocrit,snr_filter,relaxivity,steady_state_time,drift)

% A_make_R1maps_func - Generate concentration versus time curves for the
% tumor region and the arterial input region. The setup follows Loveless
% et. al 2011 MRM method of AIF filtering based upon SNR.
% 
% Inputs:
%  dce_path           - NIFTI file that contains the dynamic
%    dataset. Formated as volume x time points.
%  t1_aif_path        - NIFTI file that contains the T1 map of
%    the arterial input region (or the reference region in the case of
%    the reference region method). All values outside of the AIF ROI set
%    to <=0
%  t1_roi_path        - NIFTI file that contains the T1 map of
%    the DCE region, the region that will have DCE values calculated.
%    All values outside of the ROI set to <=0
%  noise_path         - NIFTI file that delineates a noise region of the
%    image (not T1 map). Used for SNR calculation. All values outside of
%    the ROI set to <=0
%  tr                 - reptition time (in ms) of dynamic scan
%  fa                 - flip angle (in degrees) of dynamic scan
%  hematocrit         - hematocrit percent (0 - 1.00) of subject
%  snr_filter         - snr required for AIF voxels, snr must exceed this
%    value averaged over all time points
%  relaxivity         - r1 relaxivity (in mmol^-1*sec^-1) of contrast agent 
%  steady_state_time  - in image number, defines the steady state period
%    before contrast injection. -1 indicates users will be prompted to
%    select it graphically
%  drift              - boolean value, perform drift correction based on
%    a rod phantom in image FOV
% 
% We assume that the T1 maps and the DCE-MRI files contain the same field of
% view.
% 
% The script uses the T1 maps to calculate the R1 changes over time for the
% dynamic dataset using the gradient echo signal equation. 
% 
% Voxels whereby the R1 calculation creates unreasonable values (e.g. complex
% values) are filtered out.
% 
% For the AIF (or reference region), a SNR filter is applied. Voxels that
% have SNR less than denoted will be filtered out. 
% 
% Concentration maps and time curves are saved. 
% 
% A matlab data .mat file is saved for further processing.
% 
% Requires:
% niftitools
% cleanAB.m
% cleanR1t.m
% 
% 
% Thomas Ng
% Caltech, Dec 2011
% Updated April 2012
% Updated Dec 2013 Sam Barnes




%% 1. Option Toggles

% Toggle to re-calculate the average values only in the regions that the
% voxel curvefit was able to give good fitting; ie. the viable regions.
viable= 0;


%% DO NOT ALTER LINES BELOW UNLESS YOU KNOW WHAT YOU ARE DOING
%% 2. a) Load the files

% Ask for file location
place = '';

[PathName1,base,ext] = fileparts(dce_path);
dynam = [base ext];

[PathName2,base,ext] = fileparts(t1_aif_path);
lv = [base ext];

[PathName3,base,ext] = fileparts(t1_roi_path);
tumor = [base ext];

[PathName4,base,ext] = fileparts(noise_path);
noise = [base ext];

% Output name
rootname = strrep(dynam,'.nii','_');

%Load DCE dataset
dynam = load_nii(fullfile(PathName1, place, dynam));

%image resolution
res  = dynam.hdr.dime.pixdim;

dynamname = dynam;
dynam = double(dynam.img);
dynamname.fileprefix



%Load TUMOR T1 map, find the voxels that encompass tumor ROI
TUMOR = load_nii(fullfile(PathName3, place, tumor));
TUMOR = double(TUMOR.img);
tumind= find(TUMOR > 0);

%Load AIF dataset
LV = load_nii(fullfile(PathName2, place, lv));
LV = double(LV.img);
lvind = find(LV > 0);

%Load noise ROI files
NOISE = load_nii(fullfile(PathName4, place, noise));
NOISE = double(NOISE.img);
noiseind = find(NOISE > 0);

% If viable toggle is on, choose the file that is generated with the
% viable tumor region.

if(viable)
    
    [VIA,PathName1,FilterIndex] = uigetfile([PathName1 '/*.nii'],'Choose Ktrans file Region file');
    VIA = load_nii(fullfile(PathName1, place, VIA));
    VIA = VIA.img;
    viaind = find(VIA ~=0);
    nonvia = setdiff(tumind, viaind);
    disp([num2str(100*numel(viaind)/numel(tumind)) '% of the tumor gives good fitting and is kept for further analysis'])
    
    tumind = viaind;
end

%Find the number of slices each volume dataset
slices = size(TUMOR, 3);

% The TR is in ms
tr = tr/1000; 

%% 3. Setup the data for next step.

warning off
% Initialize the time curve arrays
DYNAM       = [];
DYNAMLV     = [];
DYNAMNOISE  = [];
DYNAMNONVIA = [];


for i = 1:slices:size(dynam,3)
    %For each time point, we collect the dynamic information into the curve
    %arrays.
    
    currentimg       = dynam(:,:,i:i+(slices-1));
    DYNAM(end+1,:)   = currentimg(tumind);
    DYNAMLV(end+1,:) = currentimg(lvind);
    DYNAMNOISE(end+1)= std(currentimg(noiseind));
    

    % This is used for create a graphic showing the ROIs in relation to the
    % DCE MRI image.
    matchimg = currentimg;
    
    matchimg(tumind) = 300000;
    matchimg(lvind)  = 300000;
    matchimg(noiseind)=400000;
    
    if(viable)
        matchimg(nonvia) = 50000;
    end
    
    
    if(i == 1)
        nn = figure;
        
        for j = 1:slices
            subplot(2,slices, j), imagesc(currentimg(:,:,j)'), axis off
            subplot(2,slices,j+slices), imagesc(matchimg(:,:,j)'), axis off
        end
    end
    
end

% We save the ROI image as a fig, We save in the same directory as the
% dynamic file.
saveas(nn, fullfile(PathName1, [rootname 'image_ROI.fig']));

%% 4. Correcting drift

% If for some reason the MR signal drifts over time, and you have included
% a rod phantom in the FOV during the imaging, you can use the drift of the
% phantom to correct for the MR signal in the tissue of interest.

if(drift)
    for j = 1:slices
        figure(nn);
        cumatchimg = (matchimg(:,:,j));
        
        subplot(2, slices,j), title('First left click in the rod region and then noise region, right click if no rod (twice)');
        
        [x y button] = ginput(2);
        x = round(x);
        y = round(y);
        
        if(button(1) > 1)
            
            OUT = [];
        else
        OUT = findRod(dynam(:,:,j), [x(1) y(1)],[x(2) y(2)], []);
          
        title('');
        cumatchimg(OUT(:,1), OUT(:,2)) = 600000;
        
        subplot(2, slices, j+slices), imagesc(cumatchimg'), axis off
        end
        
        ROD{j}.OUT = OUT;
	end
    
    % Now we Drift correct the image
    originalimg = dynam(:,:,1:1+(slices-1));
    
    DYNAM       = [];
    DYNAMLV     = [];
    DYNAMNOISE  = [];
    DYNAMNONVIA = [];
    counter = 1;
    for i = 1:slices:size(dynam,3)
        currentimg       = dynam(:,:,i:i+(slices-1));
        for j = 1:slices 
            originalimgj = originalimg(:,:,j);
            currentimgj  = currentimg(:,:,j);
            
            % rod mean
            OUT = ROD{j}.OUT;
            
            if(~isempty(OUT))
                ind = sub2ind(size(originalimgj), OUT(:,1), OUT(:,2));
                scalefactor = mean(originalimgj(ind))/mean(currentimgj(ind));
                
                PRECORRECTED(counter,j) = mean(currentimgj(:));
                
                DRIFT(counter,j) = mean(currentimgj(ind));
                
                currentimgj = currentimgj.*scalefactor;
                
                CORRECTED(counter,j) = mean(currentimgj(:));
                
                currentimg(:,:,j) = currentimgj;
            end
        end
        
        DYNAM(end+1,:)   = currentimg(tumind);
        DYNAMLV(end+1,:) = currentimg(lvind);
        DYNAMNOISE(end+1)= std(currentimg(noiseind));
        
        if(viable)
            DYNAMNONVIA(end+1,:) = currentimg(nonvia);
        end
        
        counter = counter+1;
    end
    
    % Plot the drift for each slice
    figure,
    
    for j = 1:slices
        OUT = ROD{j}.OUT;
        
        if(~isempty(OUT))
            
            subplot(ceil(sqrt(slices)), ceil(sqrt(slices)), j), hold on,plot(DRIFT(:,j)', 'r.'),  plot(CORRECTED(:,j), 'gx'), plot(PRECORRECTED(:,j), 'b.'), hold off
            
            title(['Slice: ' num2str(j)]);
        end
    end
end


%% 5. Manually select injection point, if required
if(steady_state_time == -1)
    figure, hold on,plot(mean(DYNAM, 2), 'r'), plot(mean(DYNAMLV,2), 'b'),
    
    for i = 1:2
        title(['Select timepoint ' num2str(i) ...
			' before injection. (i.e. Select an interval (2 points) before contrast injection to define stead state)'])
        [steady_state_time(i) y] = ginput(1);
	end
else
	%No zero index in matlab
	steady_state_time(2) = steady_state_time;
	steady_state_time(1) = 1;
end

% Averaged timecurves for AIF and Tumor
RawTUM = mean(DYNAM,2);
RawLV  = mean(DYNAMLV,2);

if(viable)
    RawNONVIA = mean(DYNAMNONVIA,2);
end


%% 6.  Filter out AIF points if SNR too low

% Save the voxel IDs of the AIF pre filtering.
oldvind = lvind;
prefilter_number_aif_voxels = size(lvind,1);

snrfilter = 0;
%Dynam(time,voxels)
DYNAMLV_time_average = mean(DYNAMLV,1);
DYNAMNOISE_time_average = mean(DYNAMNOISE);
voxelSNR = DYNAMLV_time_average./DYNAMNOISE_time_average;

ind = find(voxelSNR < snr_filter);

lvind(ind) = [];
DYNAMLV(:,ind) = [];
voxelSNR_filtered = voxelSNR;
voxelSNR_filtered(ind) = [];

snrfilter = snrfilter + numel(ind);  
% for i = 1:size(DYNAM,1)
%     currentSNR = DYNAMLV(i,:)./DYNAMNOISE(i);
%     
%     ind = find(currentSNR < snr_filter);
%     
%     lvind(ind) = [];
%     DYNAMLV(:,ind) = [];
%     
%     numel(ind);
%     snrfilter = snrfilter + numel(ind);  
% end

disp(['AIF SNR filter requires average SNR > ' num2str(snr_filter)]);
disp(['AIF SNR filter has removed: ' num2str(snrfilter) ' of ' num2str(prefilter_number_aif_voxels) ' AIF voxels.']);
disp(['After filter average AIF SNR (all voxels, all time points) = ' num2str(mean(voxelSNR_filtered))]);
if snrfilter>=prefilter_number_aif_voxels
	error('Error SNR filter removed all AIF points, lower SNR filter requirement');
end

%% 7. Convert AIF/ Reference region to R1
T1LV     = LV(lvind);
T1       = T1LV;
Stotallv = DYNAMLV;


% Sss is the steady state Signal before injection
Sss      = mean(Stotallv(round(steady_state_time(1)):round(steady_state_time(2)),:));
Sstar    = ((1-exp(-tr./T1))./(1-cosd(fa).*exp(-tr./T1)));
Stlv     = Stotallv;%(inj(2):end,:);

for j = 1:size(Stlv,2)
    
    A(:,j) = 1-cosd(fa).*Sstar(j).*Stlv(:,j)./Sss(j);
    B(:,j)  = 1-Sstar(j).*Stlv(:,j)./Sss(j);
    
end

AB = A./B;
% return
% AB should not be less than 0. We interpolate the timeseries to clean
% this. Threshold is 0.5;
% up.
[AB T1LV lvind BADspacelv GOODspacelv] = cleanAB(AB, T1LV,lvind, 'AIF', min(numel(T1LV)), 0.5);

R1tLV = double((1/tr).*log(AB));

% R1 should be real. We interpolate the timeseries to clean
% this. Threshold is 0.5;
% up.
[R1tLV T1LV lvind BADspacelv GOODspacelv] = cleanR1t(R1tLV, T1LV,lvind, 'AIF', min(numel(T1LV)), 0.5);

for j = 1:numel(T1LV)
    
    % Scale Sss values to T1 values
    
    ScaleFactorlv = (1/T1LV(j)) - mean(R1tLV(round(steady_state_time(1)):round(steady_state_time(2)), j));
    
    R1tLV(:,j) = R1tLV(:,j) + ScaleFactorlv;
end

%% 8. Convert Tumor ROI to R1

T1TUM     = TUMOR(tumind);
T1        = T1TUM;
Stotaltum = DYNAM;

Ssstum = mean(Stotaltum(round(steady_state_time(1)):round(steady_state_time(2)),:));
Sstar    = ((1-exp(-tr./T1))./(1-cosd(fa).*exp(-tr./T1)));
Sttum = Stotaltum;

for j = 1:numel(T1)
    
    A(:,j) = 1-cosd(fa).*Sstar(j).*Sttum(:,j)./Ssstum(j);
    B(:,j)  = 1-Sstar(j).*Sttum(:,j)./Ssstum(j);
    
end

AB = A./B;
[AB T1TUM tumind BADspace GOODspace] = cleanAB(AB, T1TUM,tumind, 'Tumor', min(numel(T1TUM)), 0.7);

R1tTOI = double((1/tr).*log(AB));

[R1tTOI T1TUM tumind BADspace GOODspace] = cleanR1t(R1tTOI, T1TUM,tumind, 'Tumor', min(numel(T1TUM)), 0.7);

for j = 1:numel(T1TUM)
    % Scale Sss values to T1 values
    ScaleFactortum = (1/T1TUM(j)) - mean(R1tTOI(round(steady_state_time(1)):round(steady_state_time(2)), j));
    R1tTOI(:,j) = R1tTOI(:,j) + ScaleFactortum;
end

n = figure; subplot(421),plot(mean(R1tTOI,2), 'r.'), title('R1 maps over image reps pre-filtering Tumor'), subplot(422), plot(mean(R1tLV,2), 'b.'), title('R1 maps over image reps pre-filtering AIF')

%% 9. Convert to concentrations

% f.1 AIF
%if(~iterremove)
Cp = (R1tLV-repmat((1./T1LV)', [size(R1tLV, 1) 1]))./(relaxivity*(1-hematocrit));
%end
% f.2 Tumor
Ct = (R1tTOI-repmat((1./T1TUM)', [size(R1tTOI, 1) 1]))./relaxivity;

%% 9.1 Filter Cp more, to include voxels whose max is within the 25-75% of
%% the mean, this was used for further filtering.

% injtime = find(mean(Cp,2) == max(mean(Cp,2)));
% MAXLV = max(Cp(injtime-100:injtime+100,:), [],1);
% IQR = iqr(MAXLV);
% ind = find(abs(MAXLV - mean(MAXLV)) > IQR);
% lvind(ind) = [];
% Cp(:,ind)    = [];
% disp(['Interquartile filter had filtered out: ' num2str(numel(ind)) ' voxels.']);

%% 10. Save the Ct file as a dynamic curve

CTFILE = zeros(size(dynam));

for i = 1:size(Ct,1)
    currentCT = zeros(size(dynam,1), size(dynam,2), slices);
    currentCT(tumind) = Ct(i,:);
    CTFILE(:,:,(i-1)*slices+[1:slices]) = currentCT;
end

% Save as dynamic file, with the TUMOR ROI only.

CC = make_nii(CTFILE, res(2:4), [1 1 1]);

save_nii(CC, fullfile(PathName1, [rootname 'dynamicCt.nii']));

%% 11. Save as delta R1 values (May be useful if the Contrast agent has longer correlation time).

deltaR1LV = R1tLV-repmat(mean(R1tLV(round(steady_state_time(1)):round(steady_state_time(2)), :), 1), [size(R1tLV,1) 1]);

deltaR1TOI= R1tTOI-repmat(mean(R1tTOI(round(steady_state_time(1)):round(steady_state_time(2)), :), 1), [size(R1tTOI,1) 1]);


%% 12. Plot the time curves

figure(n), subplot(423),plot(mean(Ct,2), 'r'), title('Ct maps over image reps pre-filtering Tumor'), subplot(424), plot(mean(Cp,2), 'b'), title('Cp maps over image reps pre-filtering LV')
figure(n), subplot(425), plot(RawTUM, 'r'),

subplot(426), plot(RawLV, 'b');

subplot(428), plot(mean(deltaR1LV,2), 'r.'), title('Delta R1 blood'), subplot(427), plot(mean(deltaR1TOI,2), 'b.'), title('Delta R1 Tissue')
saveas(n,fullfile(PathName1, [rootname 'timecurves.fig']));
%% 12. Save the file for the next Step

save(fullfile(PathName1, ['A_' rootname 'R1info.mat']));
results = fullfile(PathName1, ['A_' rootname 'R1info.mat']);

disp('Finished A');































