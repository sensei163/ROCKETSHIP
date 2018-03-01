function [meanAIF,meanSignal] = AIF_manual(image_array,concentration_array, AIF_mask) 

% AIF_MANUAL is a function that computes the AIF of a user selected region
% of the image. This is done for one slice, with the remainder of the image
% being cleared to values of 0.  

% INPUTS: 

% CONCENTRATION_ARRAY: This is the image array that has been transformed
% into concentration values. 

% AIF_MASK: The user selected region of the image that contains a what
% a 'good' vessel.  This array must be the same size as the concentration
% array. 

% OUTPUT: 

% MEAN_AIF: The mean concentration time course of the AIF ROI.  

[dimx,dimy,dimz,dimt] = size(concentration_array); 

reshaped = false;
if ndims(concentration_array) == 3 && ndims(AIF_mask) == 3 
        %if so reshpae to be 4D matrix
        % But wait, there's a length difference between dimt and the time
        % length of: 
        dimt = dimz;
        dimz = 1;
        tAIF = length(AIF_mask(1,1,:)); 
        concentration_array = reshape(concentration_array,dimx,dimy,dimz,dimt);
        AIF_mask_reshape = reshape(AIF_mask, dimx,dimy,dimz,tAIF); 
        reshaped = true; 
        
elseif dimt == 1 && dimz == 1
        error('Input dynamic images not 4D data, or time dimension = 1'); 
end

AIF_indices = find(AIF_mask_reshape(:,:,:,1)>0); 

AIF_time_courses = zeros(numel(AIF_indices),dimt); 
current_time = []; 

for i = 1 : dimt 
    
 current_time = concentration_array(:,:,:,i); 
 AIF_time_courses(:,i) = current_time(AIF_indices);
 
end 

if reshaped
    
    [dimx2, dimy2, dimt2] = size(image_array);
    
    AIF_mask(:,:,1);
    image_indices = find(AIF_mask(:,:,1)>0);
    image_time_courses = zeros(numel(image_indices), dimt2);
    current_time_image = [];
    
    for i = 1 : dimt2
    
        current_time_image = image_array(:,:,i);
        image_time_courses(:,i) = current_time_image(image_indices);
    
    end

else %the following else block for the 4d case has never been tested, so it may not work
    display('4D case!')
    [dimx2, dimy2, dimz2, dimt2] = size(image_array);
    
    AIF_mask(:,:,:,1);
    image_indices = find(AIF_mask(:,:,:,1)>0);
    image_time_courses = zeros(numel(image_indices), dimt2);
    current_time_image = [];
    for i = 1 : dimt2
        current_time_image = image_array(:,:,i);
        image_time_courses(:,i) = current_time_image(image_indices);
    end
    
end
    
meanAIF = mean(AIF_time_courses,1);
meanAIF = meanAIF';     % Transpose for subsequent operations. 

meanSignal = mean(image_time_courses,1);
meanSignal = meanSignal';

%   END OF FUNCTION 



