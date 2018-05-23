function [meanAIF,meanSignal, baseline, baseline_array] = AIF_manual(image_array,concentration_array, AIF_mask, base_time_vect, base_concentration_array) 

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

if ndims(concentration_array) == 3 %3D case 
    [dimx, dimy, dimt] = size(concentration_array);
    
    %find significant AIF values 
    AIF_indices = find(AIF_mask(:,:,1) > 0);
    AIF_time_courses = zeros(numel(AIF_indices),dimt);
    current_time = [];
    
    %place these values for each slice into an array
    for i = 1:dimt
       
        current_time = concentration_array(:,:,i);
        AIF_time_courses(:,i) = current_time(AIF_indices);
    
    end
    
    %Find significant image values
    image_indices = AIF_indices;
    image_time_courses = zeros(numel(image_indices), dimt);
    current_time_image = [];
    
    %place these values into an array
    for i = 1 : dimt
    
        current_time_image = image_array(:,:,i);
        image_time_courses(:,i) = current_time_image(image_indices);
    
    end
    
    %find the concentration values that correspond to the base_line AIF  
    base_indices = AIF_indices;
    base_concentration_time_courses = zeros(numel(base_indices),numel(base_time_vect));
    
    for i = 1 : numel(base_time_vect)
        
        base_concentration_current = base_concentration_array(:,:,i);
        base_concentration_time_courses(:,i) = base_concentration_current(base_indices);
        
    end
    
else %4D case 
    
     [dimx, dimy, dimz, dimt] = size(concentration_array);
    
    %find significant AIF values 
    AIF_indices = find(AIF_mask(:,:,:,1) > 0);
    AIF_time_courses = zeros(numel(AIF_indices),dimt);
    current_time = [];
    
    %place these values for each slice into an array
    for i = 1:dimt
       
        current_time = concentration_array(:,:,:,i);
        AIF_time_courses(:,i) = current_time(AIF_indices);
    
    end
    
    %Find significant image values
    image_indices = AIF_indices;
    image_time_courses = zeros(numel(image_indices), dimt);
    current_time_image = [];
    
    %place these values into an array
    for i = 1 : dimt
    
        current_time_image = image_array(:,:,:,i);
        image_time_courses(:,i) = current_time_image(image_indices);
    
    end
    
    %find the concentration values that correspond to the base_line AIF  
    base_indices = AIF_indices;
    base_concentration_time_courses = zeros(numel(base_indices),numel(base_time_vect));
    
    for i = 1 : numel(base_time_vect)
        
        base_concentration_current = base_concentration_array(:,:,:,i);
        base_concentration_time_courses(:,i) = base_concentration_current(base_indices);
        
    end
end


meanAIF = mean(AIF_time_courses,1);
meanAIF = meanAIF';     % Transpose for subsequent operations. 

meanSignal = mean(image_time_courses,1);
meanSignal = meanSignal';

baseline_array = mean(base_concentration_time_courses)';
baseline = mean(mean(base_concentration_time_courses));

%   END OF FUNCTION 



