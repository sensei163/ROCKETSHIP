%% Clean up the R1 by interpolation

function [R1t T1 ROIind BADspace GOODspace] = cleanR1t(R1t, T1, ROIind, type, pick, threshold);

BADspace = [];
CLEANspace = 0;
timepoints = size(R1t,1)
spacepoints= size(R1t,2)

% Search each voxel to see if issues, interpolate if issues

for j = 1:spacepoints
    
    TEST = R1t(:, j);
    
    TEST = imag(TEST);
    
    ind  = find(TEST > 0);
    
    % if greater than threshold, remove
    if(ind > threshold*(timepoints))
        BADspace = [BADspace j];
    elseif(ind > 0)
        
        for k = 1:numel(ind)
            
            buffer = 3;
            
            bufferind = ind(k)-buffer:ind(k)+buffer;
            
            out = find(bufferind > (timepoints));
            bufferind(out) = [];
            
            out = find(bufferind < 1);
            bufferind(out) = [];
            
            out = find(bufferind == ind(k));
            bufferind(out) = [];
            
            CLEANspace = CLEANspace+1;
            
            
        end

        TEST(ind(k)) = interp1(bufferind, TEST(bufferind), ind(k));
        
        R1t(:, j) = TEST;
    else
    end
end

GOODspace= setdiff([1:spacepoints], BADspace);
perbad = numel(BADspace)/spacepoints;

disp([num2str(perbad) ' of total ' type ' voxels were removed from analysis.']);

disp([num2str(CLEANspace/numel(R1t)) ' of total voxels were cleaned.']);

R1t(:, BADspace) = [];
T1(BADspace)    = [];
ROIind(BADspace) = [];

% Now we pick number of voxels to be selected for actual measure 0 if all
pick = min(pick, numel(ROIind));
if(pick <= 0)
    %All chosen
else

    R1t = R1t(:,1:pick);
    
    T1 = T1(1:pick);
    ROIind = ROIind(1:pick);
end

        


        
        
        
        
        
        
        
        
        