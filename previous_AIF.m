function [meanAIF_adjusted, time_vect, concentration_array] = previous_AIF(meanAIF, meanSignal, bolus_time, time_vect, concentration_array)
%AIF_IMPORT: creates an AIF of the correct time length from an imported AIF
%Adjusted the time vector if neccessary
%   Finds the shorter array (meanAIF (imported) or time_vect) and then take elements
%   from the imported AIF 
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

