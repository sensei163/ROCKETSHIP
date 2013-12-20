function x = FXLfit_novpC(xdata, numvoxels)
  %numvoxels = xdata{1}.numvoxels;
%matlabpool open local 7

parfor i = 1:numvoxels
    
  x(i,:) =  FXLStep1AIFhelper(xdata, i)

%disp('moo')
end

%matlabpool close

