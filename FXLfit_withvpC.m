%% Helper for C_fitwithvp

function GG = FXLfit_withvpC(xdata, numvoxels)
%i = numvoxels;
%numvoxels = xdata{1}.numvoxels;
parfor i = 1:numvoxels
    
    GG(i,:) = FXLStep1AIFhelper_vp(xdata,i);
end

end
