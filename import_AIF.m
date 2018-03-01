function [meanAIF_adjusted, time_vect, concentration_array] = import_AIF(meanAIF, meanSignal, bolus_time, time_vect, concentration_array)
%AIF_IMPORT: creates an AIF of the correct time length from an imported AIF
%Adjusted the time vector if neccessary
%   Finds the shorter array (meanAIF (imported) or time_vect) and then take elements
%   from the imported AIF 
    
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




 