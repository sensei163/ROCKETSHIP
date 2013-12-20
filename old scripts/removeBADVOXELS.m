function [A B T1map ROIind] = removeBADVOXELS(A, B, T1map, ROIind, type, pick)

% This function checks the time curves A and B and remove those which are 0
% or negative to ensure reliable R1 curves are derived 

[timepoints spacepoints] = find((A <= 0) & (B>=0));

[timepoints1 spacepoints1] = find((B <=0) & (A>=0));

[timepoints2 spacepoints2] = find(abs(A) < abs(B));

BADspace = unique(union(unique([spacepoints; spacepoints1]), spacepoints2));

perbad = numel(BADspace)/size(A,2);

disp([num2str(perbad) ' of total ' type ' voxels were removed from analysis.']); 

A(:, BADspace) = [];
B(:, BADspace) = [];

T1map(BADspace) = [];
ROIind(BADspace) = [];

% Now we pick number of voxels to be selected for actual measure 0 if all
pick = min(pick, size(A,2));
if(pick <= 0)
    %All chosen
else
    
    A = A(:,1:pick);
    B = B(:,1:pick);
    
    T1map = T1map(1:pick);
    ROIind = ROIind(1:pick);
end
    

