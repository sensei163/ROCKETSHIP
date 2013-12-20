function [A B T1map ROIind R1tLV perbad Cp] = removeBADVOXELSiter(A, B, T1map, ROIind, type, pick, TR, r1, hct, inj)

% This function checks the time curves A and B and remove those which are 0
% or negative to ensure reliable R1 curves are derived

% This differs from removeBADVOXELS in that it does this for each time
% point and derive mean R1tLV = double((1/TR).*log(A./B));


timepoints = size(A,1);

for i = 1:timepoints
    
    curA = A(i,:);
    curB = B(i,:);
    curT1= T1map;
   
    
    % Find the points that need to be excluded
    
    indA = find((curA <=0) & (curB>=0));
    indB = find((curB <=0) & (curA>=0));
    indC = find(abs(curA) < abs(curB));
    
    indA = indA(:);
    indB = indB(:);
    indC = indC(:);
    
    %size(A)
    
    E = [indA;indB];
    %size(E)
    BADspace = unique(union(unique(E), indC));
    
    perbad(i) = numel(BADspace)/size(A,2);
    
    curA(BADspace) = [];
    curB(BADspace) = [];
    curT1(BADspace)= [];
    
    curR1tLV       = (double((1/TR).*log(curA./curB)));

   %save('temper.mat')
    
    
    R1tLV(i,1)     = mean(curR1tLV);
    
    curCP = (curR1tLV-repmat((1./curT1)', [size(curR1tLV, 1) 1]))./(r1*(1-hct));
    
    curCP = curCP(curCP > 0);
    
    Cp(i,1) = mean(curCP);% mean((curR1tLV-repmat((1./curT1)', [size(curR1tLV, 1) 1]))./(r1*(1-hct)))
    
   
end

% Calibrate Cp

Cpsteadstate =mean(Cp(round(inj(1)):round(inj(2))));

Cp = Cp-Cpsteadstate;
    
    
%{
[timepoints spacepoints] = find((A <= 0) & (B>=0));

[timepoints1 spacepoints1] = find((B <=0) & (A>=0));

BADspace = unique([spacepoints; spacepoints1]);

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
    
%}
