% Helper for C_fitwithvp

function GG = FXLfit_generic(xdata, number_voxels, model)

p = ProgressBar(number_voxels);

if strcmp(model,'aif_vp');
	% Preallocate for speed
	GG = zeros([number_voxels 10],'double');
	parfor i = 1:number_voxels
		GG(i,:) = FXLStep1AIFhelper_vp(xdata,i);
		p.progress;
	end;
elseif strcmp(model,'aif');
	% Preallocate for speed
	GG = zeros([number_voxels 8],'double');
	parfor i = 1:number_voxels
		GG(i,:) = FXLStep1AIFhelper(xdata,i);
		p.progress;
	end;
else
	warning(['Error, model ' model ' not yet implemented']);
	return
end

p.stop;


end
