function saved_results = A_make_R1maps_func(DYNAMIC, LV, TUMOR, NOISE, hdr, res,quant, rootname, dynampath, dynam_name, aiforRR, ... 
    tr,fa,hematocrit,snr_filter,relaxivity,steady_state_time, drift, sliceloc)

% A_make_R1maps_func - Generate concentration versus time curves for the
% tumor region and the arterial input region. The setup follows Loveless
% et. al 2011 MRM method of AIF filtering based upon SNR.
% 
% Inputs:
%  DYNAMIC            - Image matrix that contains the dynamic
%                       dataset. Formated as volume x time points.
%  LV                 - Image matrix  that contains the T1 map of
%                       the arterial input region (or the reference region
%                       in the case of the reference region method). All
%                       values outside of the AIF ROI set to <=0
%  TUMOR              - Image matrix  that contains the T1 map of the DCE
%                       region, the region that will have DCE values
%                       calculated. All values outside of the ROI set to
%                       <=0
%  NOISE              - Image matrix  that delineates a noise region of the
%                       image (not T1 map). Used for SNR calculation. All
%                       values outside of the ROI set to <=0
%  hdr                - Image header of NIFTI or DICOM file read from DCE
%                       region
%  res                - vector with voxel size
%  quant              - boolean value indicating whether to pursue
%                       quantitative DCE vs. semi-quant values
%  rootname           - Name of output file
%  dynampath          - Name of output path
%  dynamname          - Name of final file for part D
%  aiforRR            - Toggle indicating whether the auxiliary ROI is AIF (1)
%                       or reference region (2)
%  tr                 - reptition time (in ms) of dynamic scan
%  fa                 - flip angle (in degrees) of dynamic scan
%  hematocrit         - hematocrit percent (0 - 1.00) of subject
%  snr_filter         - snr required for AIF voxels, snr must exceed this
%                       value averaged over all time points
%  relaxivity         - r1 relaxivity (in mmol^-1*sec^-1) of contrast agent
%  steady_state_time  - in image number, defines the end of the steady
%                       state period before contrast injection. -1
%                       indicates users will be prompted to select it
%                       graphically
%  drift              - boolean value, perform drift correction based on
%                       a rod phantom in image FOV
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
% A_make_R1maps_func.m
% DataHash.m
% cleanAB.m
% cleanR1t.m
% findRod.m
% niftitools
% 
% Thomas Ng
% Caltech, Dec 2011
% Updated April 2012
% Updated Dec 2013 Sam Barnes
% Updated Jan 2014 Thomas Ng




%% 1. Option Toggles

% Toggle to re-calculate the average values only in the regions that the
% voxel curvefit was able to give good fitting; ie. the viable regions.
viable= 0;


%% DO NOT ALTER LINES BELOW UNLESS YOU KNOW WHAT YOU ARE DOING
% Log input results
%[log_path,base,~] = fileparts(dce_path);
log_path = dynampath;
log_path = fullfile(log_path, ['A_' rootname 'R1info.log']);
if exist(log_path, 'file')==2
  delete(log_path);
end
diary(log_path);
fprintf('************** User Input **************\n\n');
% disp('User selected dce file: ');
% fprintf('%s\n\n',dce_path);
% disp('User selected T1 of AIF: ');
% fprintf('%s\n\n',t1_aif_path);
% disp('User selected T1 of region to process: ');
% fprintf('%s\n\n',t1_roi_path);
% disp('User selected noise region in image (not T1 map): ');
% fprintf('%s\n\n',noise_path);
disp('User selected TR (ms): ');
disp(tr);
disp('User selected FA (degrees): ');
disp(fa);
disp('User selected hematocit (0 to 1.0): ');
disp(hematocrit);
disp('User selected SNR threshold for AIF: ');
disp(snr_filter);
disp('User selected contrast agent R1 relaxivity (/mmol/sec): ');
disp(relaxivity);
disp('User selected end of steady state time (image number): ');
disp(steady_state_time);
disp('User selected drift correction: ');
disp(drift);
fprintf('************** End User Input **************\n\n\n');

disp('Starting Part A Processing')
disp(datestr(now))
disp(' ');
tic

%% 2. a) Load the files
% Ask for file location
place = '';

[PathName1,~,~] = fileparts(log_path);
% [PathName1,base,ext] = fileparts(dce_path);
% dynam = [base ext];
% 
% [PathName2,base,ext] = fileparts(t1_aif_path);
% lv = [base ext];
% 
% [PathName3,base,ext] = fileparts(t1_roi_path);
% tumor = [base ext];
% 
% [PathName4,base,ext] = fileparts(noise_path);
% noise = [base ext];
% 
% % Output name
% rootname = strrep(dynam,'.nii','_');

% %Load DCE dataset
% dynam = load_nii(fullfile(PathName1, place, dynam));
% 
% %image resolution
% res  = dynam.hdr.dime.pixdim;
% save for part D
%dynam.fileprefix;
dynam = DYNAMIC; %double(dynam.img);


%Load TUMOR T1 map, find the voxels that encompass tumor ROI
%convert T1 from ms to sec
% TUMOR = load_nii(fullfile(PathName3, place, tumor));
% TUMOR = 1/1000.*double(TUMOR.img);
tumind= find(TUMOR > 0);

testt1 = mean(TUMOR(tumind));
if testt1 > 50
    disp('T1 maps of data likely in ms, converting to s...');
    TUMOR = (1/1000).*double(TUMOR);
else
    disp('T1 maps of data likely in sec, no conversion');
end



%Load AIF dataset, convert T1 from ms to sec
% LV = load_nii(fullfile(PathName2, place, lv));
% LV = 1/1000.*double(LV.img);
lvind = find(LV > 0);

testt1 = mean(LV(lvind));
if testt1 > 50
    disp('T1 maps of AIF likely in ms, converting to s...');
    LV = (1/1000).*double(LV);
else
    disp('T1 maps of AIF likely in sec, no conversion');
end

%Load noise ROI files
% NOISE = load_nii(fullfile(PathName4, place, noise));
% NOISE = double(NOISE.img);
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
    DYNAMNOISE(end+1)= std(single(currentimg(noiseind)));
    

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
        
        [x, y, button] = ginput(2);
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
    %Gets all slices from the second time point
    originalimg = dynam(:,:,(1+slices):(1+slices)+(slices-1));
    
    DYNAM       = [];
    DYNAMLV     = [];
    DYNAMNOISE  = [];
    DYNAMNONVIA = [];
%     counter = 1;
%     for i = 1:slices:size(dynam,3)
%         currentimg       = dynam(:,:,i:i+(slices-1));
%         for j = 1:slices 
%             originalimgj = originalimg(:,:,j);
%             currentimgj  = currentimg(:,:,j);
%             
%             % rod mean
%             OUT = ROD{j}.OUT;
%             
%             if(~isempty(OUT))
%                 ind = sub2ind(size(originalimgj), OUT(:,1), OUT(:,2));
%                 scalefactor = mean(originalimgj(ind))/mean(currentimgj(ind));
%                 
%                 PRECORRECTED(counter,j) = mean(currentimgj(:));
%                 
%                 DRIFT(counter,j) = mean(currentimgj(ind));
%                 
%                 currentimgj = currentimgj.*scalefactor;
%                 
%                 CORRECTED(counter,j) = mean(currentimgj(:));
%                 
%                 currentimg(:,:,j) = currentimgj;
%             end
%         end
%         
%         DYNAM(end+1,:)   = currentimg(tumind);
%         DYNAMLV(end+1,:) = currentimg(lvind);
%         DYNAMNOISE(end+1)= std(currentimg(noiseind));
%         
%         if(viable)
%             DYNAMNONVIA(end+1,:) = currentimg(nonvia);
%         end
%         
%         counter = counter+1;
%     end

    scalefactor = ones(size(dynam,3)/slices,1);
    scale_fit = cell(1,slices);
    % Slice loop
    for j = 1:slices 
        originalimgj = originalimg(:,:,j);
        
        % Time loop
        time_index = 1;
        for i = 1:slices:size(dynam,3)
            currentimg       = dynam(:,:,i:i+(slices-1));
            currentimgj  = currentimg(:,:,j);

            % rod mean
            OUT = ROD{j}.OUT;

            if(~isempty(OUT))
                ind = sub2ind(size(originalimgj), OUT(:,1), OUT(:,2));
                scalefactor(time_index) = mean(originalimgj(ind))/mean(currentimgj(ind));
            end
            time_index = time_index+1;
        end
        scale_fit{j} = fit((1:numel(scalefactor))',scalefactor,'poly3');
    end

        
    % Time loop
    time_index = 1;
    for i = 1:slices:size(dynam,3)
        currentimg       = dynam(:,:,i:i+(slices-1));
        % Slice loop
        for j = 1:slices 
            originalimgj = originalimg(:,:,j);
            currentimgj  = currentimg(:,:,j);
            
            % rod mean
            OUT = ROD{j}.OUT;
            
            if(~isempty(OUT))
                DRIFT(time_index,j) = mean(currentimgj(ind));
                
                PRECORRECTED(time_index,j) = mean(currentimgj(:));
                
                currentimgj = currentimgj.*scale_fit{j}(time_index);
                
                CORRECTED(time_index,j) = mean(currentimgj(:));
                
                currentimg(:,:,j) = currentimgj;
            end
        end
        
        DYNAM(end+1,:)   = currentimg(tumind);
        DYNAMLV(end+1,:) = currentimg(lvind);
        DYNAMNOISE(end+1)= std(single(currentimg(noiseind)));
        
        if(viable)
            DYNAMNONVIA(end+1,:) = currentimg(nonvia);
        end
        
        time_index = time_index+1;
    end
    
    
    % Plot the drift for each slice
    figure,
    
    for j = 1:slices
        OUT = ROD{j}.OUT;
        
        if(~isempty(OUT))
            
            subplot(ceil(sqrt(slices)), ceil(sqrt(slices)), j)
            hold on
                plot(DRIFT(:,j)', 'r.')
                plot(1./scale_fit{j}(1:time_index).*DRIFT(2,j),'k')
                plot(CORRECTED(:,j), 'gx')
                plot(PRECORRECTED(:,j), 'b.')
            hold off
            
            title(['Slice: ' num2str(j)]);
        end
    end
end


%% 5. Manually select injection point, if required
if(steady_state_time == -1)
    figure, hold on;
    plot(mean(DYNAM, 2), 'r');
    plot(mean(DYNAMLV,2), 'b');
    ylabel('a.u.');
    xlabel('image number');
    hold off;
    
    
    for i = 1:2
        title(['Select timepoint ' num2str(i) ...
            ' before injection. (i.e. Select an interval (2 points) before contrast injection to define stead state)'])
        [steady_state_time(i), ~] = ginput(1);
    end
    steady_state_time(2) = round(steady_state_time(2));
    steady_state_time(1) = round(steady_state_time(1));
    if steady_state_time(1)<1
        %No zero index in matlab
        steady_state_time(1) = 1;
    end
else
    %No zero index in matlab
    steady_state_time(2) = steady_state_time;
    steady_state_time(1) = 1;
end
disp(['Steady state time selected from image ' num2str(steady_state_time(1)) ...
    ' to image ' num2str(steady_state_time(2)) ' ']);

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

n = figure; 
subplot(421),plot(mean(R1tTOI,2), 'r.'), title('R1 maps pre-filtering ROI'), ylabel('sec^-1')
subplot(422), plot(mean(R1tLV,2), 'b.'), title('R1 maps pre-filtering AIF'), ylabel('sec^-1')

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

CC = make_nii(CTFILE, res, [1 1 1]);

save_nii(CC, fullfile(PathName1, [rootname 'dynamicCt.nii']));

%% 11. Save as delta R1 values (May be useful if the Contrast agent has longer correlation time).

deltaR1LV = R1tLV-repmat(mean(R1tLV(round(steady_state_time(1)):round(steady_state_time(2)), :), 1), [size(R1tLV,1) 1]);

deltaR1TOI= R1tTOI-repmat(mean(R1tTOI(round(steady_state_time(1)):round(steady_state_time(2)), :), 1), [size(R1tTOI,1) 1]);


%% 12. Plot the time curves

figure(n), 
subplot(423),plot(mean(Ct,2), 'r'), title('Ct maps pre-filtering ROI'), ylabel('mmol')
subplot(424), plot(mean(Cp,2), 'b'), title('Cp maps pre-filtering AIF'), ylabel('mmol')
subplot(425), plot(RawTUM, 'r'), title('T1-weighted ROI'), ylabel('a.u.')
subplot(426), plot(RawLV, 'b'), title('T1-weighted AIF'), ylabel('a.u.')
subplot(427), plot(mean(deltaR1TOI,2), 'b.'), title('Delta R1 ROI'), ylabel('sec^-1')
subplot(428), plot(mean(deltaR1LV,2), 'r.'), title('Delta R1 AIF') , ylabel('sec^-1')
saveas(n,fullfile(PathName1, [rootname 'timecurves.fig']));

%% 13. Setup output structure

Adata.CTFILE = CTFILE;
Adata.Cp     = Cp;
Adata.Ct     = Ct;
Adata.DYNAM  = DYNAM;
Adata.DYNAMIC= DYNAMIC;
Adata.DYNAMLV= DYNAMLV;
Adata.DYNAMLV_time_average = DYNAMLV_time_average;
Adata.DYNAMNOISE = DYNAMNOISE;
Adata.DYNAMNOISE_time_average = DYNAMNOISE_time_average;
Adata.DYNAMNONVIA = DYNAMNONVIA;
Adata.GOODspace = GOODspace;
Adata.GOOODspacelv = GOODspacelv;
Adata.LV = LV;
Adata.NOISE = NOISE;
Adata.PathName1 = PathName1;
Adata.R1tLV = R1tLV;
Adata.R1tTOI= R1tTOI;
Adata.RawLV = RawLV;
Adata.RawTUM = RawTUM;
Adata.ScaleFactorlv = ScaleFactorlv;
Adata.ScaleFactortum= ScaleFactortum;
Adata.Sss = Sss;
Adata.Ssstum = Ssstum;
Adata.Sstar  = Sstar;
Adata.Stlv   = Stlv;
%Adata.Stotaltum = Stotaltum;
Adata.Sttum  = Sttum;
Adata.T1 = T1;
Adata.T1LV = T1LV;
Adata.T1TUM=T1TUM;
Adata.TUMOR= TUMOR;
Adata.aiforRR = aiforRR;
Adata.currentCT=currentCT;
Adata.deltaR1LV = deltaR1LV;
Adata.deltaR1TOI= deltaR1TOI;
Adata.dynam_name= dynam_name;
Adata.dynampath = dynampath;
Adata.hdr       = hdr;
Adata.hematocrit=hematocrit;
Adata.log_path  = log_path;
Adata.matchimg  = matchimg;
Adata.noiseind  = noiseind;
Adata.oldvind   = oldvind;
Adata.quant     = quant;
Adata.relaxivity= relaxivity;
Adata.res       = res;
Adata.rootname  = rootname;
Adata.sliceloc  = sliceloc;
Adata.steady_state_time = steady_state_time;
Adata.tumind    = tumind;
Adata.voxelSNR  = voxelSNR;
Adata.voxelSNR_filtered = voxelSNR_filtered;

%% 14. Save the file for the next Step


saved_results = fullfile(PathName1, ['A_' rootname 'R1info.mat']);
save(saved_results, 'Adata');
Opt.Input = 'file';
mat_md5 = DataHash(saved_results, Opt);
disp(' ')
disp('MAT results saved to: ')
disp(saved_results)
disp(['File MD5 hash: ' mat_md5])

disp(' ');
disp('Finished A');
disp(datestr(now))
toc
diary off;





























