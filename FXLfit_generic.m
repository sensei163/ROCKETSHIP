% Helper for C_fitwithvp

function GG = FXLfit_generic(xdata, number_voxels, model)

p = ProgressBar(number_voxels);

if strcmp(model,'aif_vp');
	% Preallocate for speed
	GG = zeros([number_voxels 10],'double');
	% Slice out needed variables for speed
	Ct_data = xdata{1}.Ct;
	Cp_data = xdata{1}.Cp;
	timer_data = xdata{1}.timer;
	parfor i = 1:number_voxels
		GG(i,:) = FXLStep1AIFhelper_vp(Ct_data(:,i),Cp_data,timer_data);
		p.progress;
	end;
elseif strcmp(model,'aif');
	% Preallocate for speed
	GG = zeros([number_voxels 8],'double');
	% Slice out needed variables for speed
	Ct_data = xdata{1}.Ct;
	Cp_data = xdata{1}.Cp;
	timer_data = xdata{1}.timer;
	parfor i = 1:number_voxels
		GG(i,:) = FXLStep1AIFhelper(Ct_data(:,i),Cp_data,timer_data);
		p.progress;
	end;
else
	warning(['Error, model ' model ' not yet implemented']);
	return
end

p.stop;


end
