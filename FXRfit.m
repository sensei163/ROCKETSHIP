function GG = FXRfit(xdata, numvoxels)
 % numvoxels = xdata{1}.numvoxels;
%matlabpool open local 7
		  %


parfor i = 1:numvoxels
    
  GG(i,:) =  FXRAIFhelper(xdata, i);

end


 end
%matlabpool close

