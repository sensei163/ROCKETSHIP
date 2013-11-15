% Helper for C_fitwithvp

function GG = FXLfit_withvpC(xdata, number_voxels)
%i = numvoxels;
%numvoxels = xdata{1}.numvoxels;

% Preallocate for speed
GG = zeros([number_voxels 4],'double');

p = ProgressBar(number_voxels);
parfor i = 1:number_voxels
    GG(i,:) = FXLStep1AIFhelper(xdata,i);
	p.progress;
end
p.stop;


end
