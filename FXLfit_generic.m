% Helper for C_fitwithvp

function GG = FXLfit_generic(xdata, number_voxels, model, smooth_model)

% Preallocate for speed
GG = zeros([number_voxels 4],'double');

p = ProgressBar(number_voxels);

if strcmp(model,'aif_vp');
	parfor i = 1:number_voxels
		GG(i,:) = FXLStep1AIFhelper_vp(xdata,i, smooth_model);
		p.progress;
	end;
elseif strcmp(model,'aif');
	parfor i = 1:number_voxels
		GG(i,:) = FXLStep1AIFhelper(xdata,i, smooth_model);
		p.progress;
	end;
else
	warning(['Error, model ' model ' not yet implemented']);
	return
end

p.stop;


end
