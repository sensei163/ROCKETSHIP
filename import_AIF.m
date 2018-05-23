function [meanAIF_adjusted, time_vect, concentration_array, meanSignal] = import_AIF(meanAIF, bolus_time, time_vect, concentration_array, r2_star, TE)
%AIF_IMPORT: creates an AIF of the correct time length from an imported AIF
%Adjusted the time vector if neccessary
%   Finds the shorter array (meanAIF (imported) or time_vect) and then take elements
%   from the imported AIF 

%this code loads the file
%{
[filename, pathname, filterspec] = uigetfile('*.nii','Nifti Files (*.nii)'); 

if isequal(filename,0)
    %disp('User selected Cancel')
    fullpath = '';
else
    %disp(['User selected ', fullfile(pathname, filename)])
    
    % Combine path and filename together
    fullpath = strcat(pathname,filename);
    path_out = fullpath; 
    
    if ischar(fullpath)
        fullpath = {fullpath};
    end
    % if isempty(list)
    %     list = {''};
    % end
    
    %filename = filename';
    fullpath = fullpath';
    
    if numel(fullpath) > 1
        visualpath = [fullpath{1} ' -> ' num2str(numel(fullpath)) ' files'];
    else
        visualpath = fullpath{1};
    end
    
    set(handles.import_aif,'String',visualpath);
    handles.import_stored = path_out; 
   
end
%}
    
    %adjust mean AIF so it is from the injection time onwards
    meanAIF = meanAIF(bolus_time:end);
    if numel(time_vect) <= numel(meanAIF) %truncate meanAIF array if longer and copy into AIF_adjusted
        meanAIF_adjusted = zeros(numel(time_vect),1);
        for i = 1:numel(time_vect)
            meanAIF_adjusted(i) = meanAIF(i);
        end
    elseif numel(time_vect) > numel(meanAIF) %truncate time array and concentration array if longer than AIF
        display('Image Data Truncated to Match AIF')
        meanAIF_adjusted = zeros(numel(meanAIF),1);
        for i = 1:numel(meanAIF)
            meanAIF_adjusted(i) = meanAIF(i);
        end
        
        %truncate the time vector
        time_vect = time_vect(1:numel(meanAIF));
        
        %truncate the concentration array
        if ndims(concentration_array) == 4
            concentration_array(:,:,:,(numel(meanAIF)+1):end) = [];
        else
            concentration_array(:,:,(numel(meanAIF)+1):end) = [];           
        end
    end
    
    e = exp(1);
    
    meanSignal = e.^(-r2_star*meanAIF*TE); %for ablsolute meanSignal: mS = S0 * e.^(-r2_star*meanAIF*TE)
    %as S0 is not known, this meanSignal will only give 
    %relative meanSignal curve 
    
end

%{ 
 
AIF Input File Formate: 
    
    bolus_time is a slice number or index!!

    bolus_time = contrast injection index - start index of AIF 
    
    ex: if the contrast is injected at the same time as the AIF then the 
    bolus_time is zero 

        if AIF starts at time zero then bolus_time = start index of the AIF
        
        if the AIF starts at a nonzero index then use the formula 
        bolus_time = contrast injection index - start index of AIF to get
        the appropriate bolus_time value 
%} 




 