function [ concentration_array,base_concentration_array, time_vect, base_time_vect, whole_time_vect, bolus_time] = DSC_signal2concentration(image_array,TE,TR,r2_star,species,path,noise_type,noise_array)
% FUNCTION PURPOSE: To convert signal values to units of concentration, mM.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%INPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% IMAGE_ARRAY: The image arranged into an array of signal values. Must be a
% .nii file.
%
% TE: The echo-time used in the DSC scan.
%
% TR: The time resolution of the DSC MRI scan.
%
% r2_star:  The relaxivity value of the concentrast agent. This is
% dependent on both species and magnet used.
% 
% SPECIES: Which species is the data obtain from.
%
% PATH: The file path where the image array 'lives'.
%
% NOISE TYPE: The noise is handled two different ways. 0 indicates that the
% noise region is automatically taken from the lower right corner. 1
% indicates that the user has provided a noise region file into the GUI
% input.
%
% NOISE ARRAY: (optional) The noise region file; must be a .nii file that
% is the same size as the image arra.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%OUTPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONCENTRATION_ARRAY: An array of concentration values. This array is used
% as an input for subsequent functios.
%
% TIME_VECT a vector of time values. This will be used in subsequent
% calculations.
% 
% SAVED FILE: A nifti image, ending in nii_conc.nii, is saved under the
% same path as the image file


%%%%%%%%%%%%%%%%******BEGIN FUNCTION BODY******%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Doing some array housekeeping: 

% We want our array in double format: 
image_array = double(image_array); 

%Calculate size of dimensions for subsequent calculations:
if ndims(image_array) == 4
[dimx,dimy,dimz,dimt] = size(image_array);  
else 
[dimx,dimy,dimt] = size(image_array); 
dimz = 1; 
end

% A vector to be filled with mean signal values:
mean_signal = zeros(dimt,1);

%% Find the bolus injection time: 

% Going to compute the mean and look for a 'dip'. This will be done an a
% central slice. However, if we are using an image that consists of a single
% slice, this convention is neglected.

% Checking the image array and computing the mean value of select voxels for each time step. 
% This method should allow us to find the bolus time and select the proper
% time window for the subsequent CBF and CBV computations. 

if ndims(image_array) == 4 
    
    central_slice = round(dimz/2); 
    
    % does our image contain negative values: 
    if mean(mean(image_array(:,:,central_slice,1))) < 0  % If yes, then set the minimum value to 0 (make it unsigned) 
        image_array = image_array + 32768; 
    end
    % We need to find voxels scan for the bolus time. Since DSC uses a
    % changes in T2* in order to generate its signal values, we are going to look areas
    % that appear bright with respect to the background. 
    
    % We are going to take the background to to be the mean, baseline
    % intensity over the frames 2-8. 
    
    if noise_type == 0 
        background_signal = mean(mean(mean(image_array(end-4:end,end-4:end,central_slice,2:8)))); 
    elseif noise_type == 1 
        noise_voxels = find(squeeze(noise_array(:,:,1))>0); 
        noise_time_courses = zeros(numel(noise_voxels),7); 
        current_time = []; 
        for i = 1 : 7
            current_time = squeeze(noise_array(:,:,central_slice,i+1)); 
            noise_time_courses(:,i) = current_time(noise_voxels); 
        end
        background_signal = mean(mean(noise_time_courses)); 
    end
        
     if strcmp(species,'mouse')
         [rows, cols] = find(image_array(:,:,central_slice,2) > background_signal * 3); 
     elseif strcmp(species,'human') 
         [rows, cols] = find(image_array(:,:,central_slice,2) > 3000); 
     end
         
 
    % Now we check only the elements dicated by the output of FIND above: 
    check_array = zeros(length(rows),1); 
    for i = 1 : dimt
        for j = 1 : length(rows) 
                check_array(j,1) = image_array(rows(j),cols(j),central_slice,i); 
        end
        mean_signal(i) = mean(check_array); 
    end
    
else % We are working with a slice
    
    if mean(mean(image_array(:,:,1))) < 0
        image_array = image_array + 32768; 
    end
    
    if noise_type == 0
       background_signal = mean(mean(mean(image_array(end-4:end,end-4:end,2:8)))); 
    elseif noise_type == 1 
                noise_voxels = find(squeeze(noise_array(:,:,1))>0); 
        noise_time_courses = zeros(numel(noise_voxels),7); 
        current_time = []; 
        for i = 1 : 7
            current_time = noise_array(:,:,i+1); 
            noise_time_courses(:,i) = squeeze(current_time(noise_voxels)); 
        end
        background_signal = mean(mean(noise_time_courses)); 
    end
        

     if strcmp(species,'mouse')
         [rows, cols] = find(image_array(:,:,2) > background_signal * 3); 
     elseif strcmp(species,'human') 
         [rows, cols] = find(image_array(:,:,2) > 3000); 
     end
         
    % Now we check only the elements dicated by the output of FIND above: 
    check_array = zeros(length(rows),1); 
    for i = 1 : dimt
        for j = 1 : length(rows) 
                check_array(j,1) = image_array(rows(j),cols(j),i); 
        end
        mean_signal(i) = mean(check_array); 
    end
end

% Now we need to find the bolus time. Looking for a 'significant' dip! 
stdev_sigs = zeros(dimt-2,1); 
[~,minsignal_index] = min(mean_signal);
for i = 3 : dimt 
    stdev_sigs(-2+i,1) = std(mean_signal(3:i)); 
    if (i > 6) && ( (mean_signal(i)) < (mean(mean_signal(3:i-1)) - 4*mean(stdev_sigs(2:-2+i-1,1))))...
            && ((minsignal_index - i) <= 4) %added 5/21/18 to prevent bolus_time selection during the baseline phase
        bolus_time = i-1;
        break
    end 
end 

%% Calculate baseline signal before injection and compute concentration: 

% Initialize an array of concentration values: 
if ndims(image_array) == 4
    concentration_array = zeros(dimx,dimy,dimz,dimt-bolus_time); 
    % We need an array for the baseline values, too: 
    base_signal_array = zeros(dimx,dimy,dimz); 
else 
    concentration_array = zeros(dimx,dimy,dimt-bolus_time); 
    base_signal_array = zeros(dimx,dimy); 
end

    
% Now fill the arrays and compute: 

%calculate a concentrtion array from the bolus injection time + 1 until the end
%of the scan 
for t = 1 : dimt
    for k = 1 : dimz 
        for j = 1 : dimy 
            for i = 1 : dimx 
                
                if ndims(image_array) == 4
                    if strcmp(species,'human')
                            if mean(image_array(i,j,k,:)) > 3000 && mean(image_array(i,j,k,:)) < 10000  % If the value are not in this range. Put them to zero. 
                            % The baseline signal value: 
                            base_signal_array(i,j,k) = mean(image_array(i,j,k,2:bolus_time-1)); 
                            % Compute the concentration in a voxel-wise manner: 
                            %concentration_array(i,j,k,t) = (-1/TE)*log(image_array(i,j,k,t+bolus_time-1)/base_signal_array(i,j,k)) ./ r2_star; 
                            concentration_array(i,j,k,t) = (-1/TE)*log(image_array(i,j,k,t)/base_signal_array(i,j,k)) ./ r2_star;
                            end
                    elseif strcmp(species,'mouse')
                          if mean(image_array(i,j,k,:)) > background_signal * 3 && mean(image_array(i,j,k,:)) < 50000 % If the value are not in this range. Put them to zero. 
                            % The baseline signal value: 
                            base_signal_array(i,j,k) = mean(image_array(i,j,k,2:bolus_time-1)); 
                            % Compute the concentration in a voxel-wise manner: 
                            %concentration_array(i,j,k,t) = (-1/TE)*log(image_array(i,j,k,t+bolus_time-1)/base_signal_array(i,j,k)) ./ r2_star; 
                            concentration_array(i,j,k,t) = (-1/TE)*log(image_array(i,j,k,t)/base_signal_array(i,j,k)) ./ r2_star;
                          end
                    end
                        
                        
                else 
                    if strcmp(species,'human')
                        if mean(image_array(i,j,:)) >  3000 && mean(image_array(i,j,:)) < 10000
                        base_signal_array(i,j) = mean(image_array(i,j,2:bolus_time-1));  
                        %concentration_array(i,j,t) = (-1/TE)*log(image_array(i,j,t+bolus_time-1)/base_signal_array(i,j)) ./ r2_star;   
                        concentration_array(i,j,t) = (-1/TE)*log(image_array(i,j,t)/base_signal_array(i,j)) ./ r2_star; 
                        end
                    elseif strcmp(species,'mouse')
                        if mean(image_array(i,j,:)) > background_signal * 3  && mean(image_array(i,j,:)) < 50000
                        base_signal_array(i,j) = mean(image_array(i,j,2:bolus_time-1));  
                        %concentration_array(i,j,t) = (-1/TE)*log(image_array(i,j,t+bolus_time-1)/base_signal_array(i,j)) ./ r2_star; 
                        concentration_array(i,j,t) = (-1/TE)*log(image_array(i,j,t)/base_signal_array(i,j)) ./ r2_star;
                        end
                    end
            
                end
            end
        end
    end
end 

%create a baseline concetration array to hold the concentration values from
%start of the scan to the bolus injection time
%base_concentration_array = zeros(dimx,dimy, bolus_time);
%{
for t = 1 : bolus_time
    for k = 1 : dimz 
        for j = 1 : dimy 
            for i = 1 : dimx 
                
                if ndims(image_array) == 4
                    if strcmp(species,'human')
                            if mean(image_array(i,j,k,:)) > 3000 && mean(image_array(i,j,k,:)) < 10000  % If the value are not in this range. Put them to zero. 
                            % The baseline signal value: 
                            base_signal_array(i,j,k) = mean(image_array(i,j,k,2:bolus_time-1)); 
                            % Compute the concentration in a voxel-wise manner: 
                            base_concentration_array(i,j,k,t) = (-1/TE)*log(image_array(i,j,k,t)/base_signal_array(i,j,k)) ./ r2_star; 
                            end
                    elseif strcmp(species,'mouse')
                          if mean(image_array(i,j,k,:)) > background_signal * 3 && mean(image_array(i,j,k,:)) < 50000 % If the value are not in this range. Put them to zero. 
                            % The baseline signal value: 
                            base_signal_array(i,j,k) = mean(image_array(i,j,k,2:bolus_time-1)); 
                            % Compute the concentration in a voxel-wise manner: 
                            base_concentration_array(i,j,k,t) = (-1/TE)*log(image_array(i,j,k,t)/base_signal_array(i,j,k)) ./ r2_star; 
                          end
                    end
                        
                        
                else 
                    if strcmp(species,'human')
                        if mean(image_array(i,j,:)) >  3000 && mean(image_array(i,j,:)) < 10000
                        base_signal_array(i,j) = mean(image_array(i,j,2:bolus_time-1));  
                        base_concentration_array(i,j,t) = (-1/TE)*log(image_array(i,j,t)/base_signal_array(i,j)) ./ r2_star;   
                        end
                    elseif strcmp(species,'mouse')
                        if mean(image_array(i,j,:)) > background_signal * 3  && mean(image_array(i,j,:)) < 50000
                        base_signal_array(i,j) = mean(image_array(i,j,2:bolus_time-1));  
                        base_concentration_array(i,j,t) = (-1/TE)*log(image_array(i,j,t)/base_signal_array(i,j)) ./ r2_star;  
                        end
                    end
            
                end
            end
        end
    end
end 
%}
% In case we have bad values: 
concentration_array(concentration_array<0) = 0; 
concentration_array(concentration_array >32767) = 0; 
concentration_array(isnan(concentration_array)) = 0; 

%extract the concentration_array and base_concentration array 
%in this verson of the program the concentration_array is expanded to include the bolus time
%this matches the equation for the calcuation of CBV in Ostergaard et al.,
%1996)
if ndims(concentration_array) == 3
    base_concentration_array = concentration_array(:,:,1:bolus_time-1);
    concentration_array = concentration_array(:,:,(bolus_time):end);
else
    base_concentration_array = concentration_array(:,:,:,1:bolus_time-1);
    concentration_array = concentration_array(:,:,:,(bolus_time):end); 
    
end
% Now make the time vector: 

time_vect = 0 : TR : (dimt-bolus_time) * TR; 
time_vect = (time_vect ./ 60)';     % transpose and put in terms of minutes.  
%time_vect_plot = 0 : TR : (dimt-1) * TR;
base_time_vect = 0: TR : (bolus_time - 2) * TR; %
base_time_vect = (base_time_vect ./ 60)'; %transpose and put in terms of minutes

whole_time_vect = 0: TR : (dimt - 1) *TR;
whole_time_vect = (whole_time_vect ./ 60)';

%Lastly, save the concentration array as nifti image file 

nii_conc = make_nii(concentration_array); 

% Lastly, we are going to save the nifti file to the same location as the
% root image.  
file_name = strcat(path,'nii_conc.nii'); 
save_nii(nii_conc, file_name); 

%%%%%%%%%%%%%%%%%%%%%******END OF FUNCTION******%%%%%%%%%%%%%%%%%%%%%%%%%%%
