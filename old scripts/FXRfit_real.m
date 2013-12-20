%% Helper for C_fit_FXR

function GG = FXRfit_real(xdata)
  numvoxels = xdata{1}.numvoxels;
%matlabpool open local 7

parfor i = 1:2;%numvoxels
    
  GG(i,:) =  FXRAIFhelper(xdata, i)

end
end

%matlabpool close

