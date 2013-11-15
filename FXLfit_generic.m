% Helper for C_fitwithvp

function GG = FXLfit_generic(xdata, number_voxels, model)

% Preallocate for speed
GG = zeros([number_voxels 4],'double');

p = ProgressBar(number_voxels);
parfor i = 1:number_voxels
	if strcmp(model,'vp');
		GG(i,:) = FXLStep1AIFhelper_vp(xdata,i);
	elseif strcmp(model,'novp');
		GG(i,:) = FXLStep1AIFhelper(xdata,i);
	else
	end
	p.progress;
end
p.stop;


end
