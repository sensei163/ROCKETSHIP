% Helper for C_fitwithvp

function GG = FXLfit_generic(xdata, number_voxels, model)

if strcmp(model,'aif_vp');
    % Get values from pref file
    prefs_str = parse_preference_file('dce_preferences.txt',0,...
        {'voxel_lower_limit_ktrans' 'voxel_upper_limit_ktrans' 'voxel_initial_value_ktrans' ...
        'voxel_lower_limit_ve' 'voxel_upper_limit_ve' 'voxel_initial_value_ve' ...
        'voxel_lower_limit_vp' 'voxel_upper_limit_vp' 'voxel_initial_value_vp' ...
        'voxel_TolFun' 'voxel_TolX' 'voxel_MaxIter' 'voxel_MaxFunEvals' 'voxel_Robust'});
    prefs.lower_limit_ktrans = str2num(prefs_str.voxel_lower_limit_ktrans);
    prefs.upper_limit_ktrans = str2num(prefs_str.voxel_upper_limit_ktrans);
    prefs.initial_value_ktrans = str2num(prefs_str.voxel_initial_value_ktrans);
    prefs.lower_limit_ve = str2num(prefs_str.voxel_lower_limit_ve);
    prefs.upper_limit_ve = str2num(prefs_str.voxel_upper_limit_ve);
    prefs.initial_value_ve = str2num(prefs_str.voxel_initial_value_ve);
    prefs.lower_limit_vp = str2num(prefs_str.voxel_lower_limit_vp);
    prefs.upper_limit_vp = str2num(prefs_str.voxel_upper_limit_vp);
    prefs.initial_value_vp = str2num(prefs_str.voxel_initial_value_vp);
    prefs.TolFun = str2num(prefs_str.voxel_TolFun);
    prefs.TolX = str2num(prefs_str.voxel_TolX);
    prefs.MaxIter = str2num(prefs_str.voxel_MaxIter);
    prefs.MaxFunEvals = str2num(prefs_str.voxel_MaxFunEvals);
    prefs.Robust = prefs_str.voxel_Robust;
    %Log values used
    fprintf('lower_limit_ktrans = %s\n',num2str(prefs.lower_limit_ktrans));
    fprintf('upper_limit_ktrans = %s\n',num2str(prefs.upper_limit_ktrans));
    fprintf('initial_value_ktrans = %s\n',num2str(prefs.initial_value_ktrans));
    fprintf('lower_limit_ve = %s\n',num2str(prefs.lower_limit_ve));
    fprintf('upper_limit_ve = %s\n',num2str(prefs.upper_limit_ve));
    fprintf('initial_value_ve = %s\n',num2str(prefs.initial_value_ve));
    fprintf('lower_limit_vp = %s\n',num2str(prefs.lower_limit_vp));
    fprintf('upper_limit_vp = %s\n',num2str(prefs.upper_limit_vp));
    fprintf('initial_value_vp = %s\n',num2str(prefs.initial_value_vp));
    fprintf('TolFun = %s\n',num2str(prefs.TolFun));
    fprintf('TolX = %s\n',num2str(prefs.TolX));
    fprintf('MaxIter = %s\n',num2str(prefs.MaxIter));
    fprintf('MaxFunEvals = %s\n',num2str(prefs.MaxFunEvals));
    fprintf('Robust = %s\n',num2str(prefs.Robust));

    % Preallocate for speed
    GG = zeros([number_voxels 10],'double');
    % Slice out needed variables for speed
    Ct_data = xdata{1}.Ct;
    Cp_data = xdata{1}.Cp;
    timer_data = xdata{1}.timer;
    %Turn off diary if on as it doesn't work with progress bar
    diary_restore = 0;
    if strcmp(get(0,'Diary'),'on')
        diary off;
        diary_restore = 1;
    end
    p = ProgressBar(number_voxels);
    parfor i = 1:number_voxels
        GG(i,:) = FXLStep1AIFhelper_vp(Ct_data(:,i),Cp_data,timer_data,prefs);
        p.progress;
    end;
    p.stop;
    if diary_restore, diary on, end;
elseif strcmp(model,'aif');
    % Get values from pref file
    prefs_str = parse_preference_file('dce_preferences.txt',0,...
        {'voxel_lower_limit_ktrans' 'voxel_upper_limit_ktrans' 'voxel_initial_value_ktrans' ...
        'voxel_lower_limit_ve' 'voxel_upper_limit_ve' 'voxel_initial_value_ve' ...
        'voxel_TolFun' 'voxel_TolX' 'voxel_MaxIter' 'voxel_MaxFunEvals' 'voxel_Robust'});
    prefs.lower_limit_ktrans = str2num(prefs_str.voxel_lower_limit_ktrans);
    prefs.upper_limit_ktrans = str2num(prefs_str.voxel_upper_limit_ktrans);
    prefs.initial_value_ktrans = str2num(prefs_str.voxel_initial_value_ktrans);
    prefs.lower_limit_ve = str2num(prefs_str.voxel_lower_limit_ve);
    prefs.upper_limit_ve = str2num(prefs_str.voxel_upper_limit_ve);
    prefs.initial_value_ve = str2num(prefs_str.voxel_initial_value_ve);
    prefs.TolFun = str2num(prefs_str.voxel_TolFun);
    prefs.TolX = str2num(prefs_str.voxel_TolX);
    prefs.MaxIter = str2num(prefs_str.voxel_MaxIter);
    prefs.MaxFunEvals = str2num(prefs_str.voxel_MaxFunEvals);
    prefs.Robust = prefs_str.voxel_Robust;
    %Log values used
    fprintf('lower_limit_ktrans = %s\n',num2str(prefs.lower_limit_ktrans));
    fprintf('upper_limit_ktrans = %s\n',num2str(prefs.upper_limit_ktrans));
    fprintf('initial_value_ktrans = %s\n',num2str(prefs.initial_value_ktrans));
    fprintf('lower_limit_ve = %s\n',num2str(prefs.lower_limit_ve));
    fprintf('upper_limit_ve = %s\n',num2str(prefs.upper_limit_ve));
    fprintf('initial_value_ve = %s\n',num2str(prefs.initial_value_ve));
    fprintf('TolFun = %s\n',num2str(prefs.TolFun));
    fprintf('TolX = %s\n',num2str(prefs.TolX));
    fprintf('MaxIter = %s\n',num2str(prefs.MaxIter));
    fprintf('MaxFunEvals = %s\n',num2str(prefs.MaxFunEvals));
    fprintf('Robust = %s\n',num2str(prefs.Robust));
    
    % Preallocate for speed
    GG = zeros([number_voxels 8],'double');
    % Slice out needed variables for speed
    Ct_data = xdata{1}.Ct;
    Cp_data = xdata{1}.Cp;
    timer_data = xdata{1}.timer;
    %Turn off diary if on as it doesn't work with progress bar
    diary_restore = 0;
    if strcmp(get(0,'Diary'),'on')
        diary off;
        diary_restore = 1;
    end
    p = ProgressBar(number_voxels);
    parfor i = 1:number_voxels
        GG(i,:) = FXLStep1AIFhelper(Ct_data(:,i),Cp_data,timer_data,prefs);
        p.progress;
    end;
    p.stop;
    if diary_restore, diary on, end;
elseif strcmp(model,'fxr');
    % Get values from pref file
    prefs_str = parse_preference_file('dce_preferences.txt',0,...
        {'voxel_lower_limit_ktrans' 'voxel_upper_limit_ktrans' 'voxel_initial_value_ktrans' ...
        'voxel_lower_limit_ve' 'voxel_upper_limit_ve' 'voxel_initial_value_ve' ...
        'voxel_lower_limit_tau' 'voxel_upper_limit_tau' 'voxel_initial_value_tau' ...
        'voxel_TolFun' 'voxel_TolX' 'voxel_MaxIter' 'voxel_MaxFunEvals' 'voxel_Robust'...
        'fxr_fw'});
    prefs.lower_limit_ktrans = str2num(prefs_str.voxel_lower_limit_ktrans);
    prefs.upper_limit_ktrans = str2num(prefs_str.voxel_upper_limit_ktrans);
    prefs.initial_value_ktrans = str2num(prefs_str.voxel_initial_value_ktrans);
    prefs.lower_limit_ve = str2num(prefs_str.voxel_lower_limit_ve);
    prefs.upper_limit_ve = str2num(prefs_str.voxel_upper_limit_ve);
    prefs.initial_value_ve = str2num(prefs_str.voxel_initial_value_ve);
    prefs.lower_limit_tau = str2num(prefs_str.voxel_lower_limit_tau);
    prefs.upper_limit_tau = str2num(prefs_str.voxel_upper_limit_tau);
    prefs.initial_value_tau = str2num(prefs_str.voxel_initial_value_tau);
    prefs.TolFun = str2num(prefs_str.voxel_TolFun);
    prefs.TolX = str2num(prefs_str.voxel_TolX);
    prefs.MaxIter = str2num(prefs_str.voxel_MaxIter);
    prefs.MaxFunEvals = str2num(prefs_str.voxel_MaxFunEvals);
    prefs.Robust = prefs_str.voxel_Robust;
    prefs.fxr_fw = str2num(prefs_str.fxr_fw);
    %Log values used
    fprintf('lower_limit_ktrans = %s\n',num2str(prefs.lower_limit_ktrans));
    fprintf('upper_limit_ktrans = %s\n',num2str(prefs.upper_limit_ktrans));
    fprintf('initial_value_ktrans = %s\n',num2str(prefs.initial_value_ktrans));
    fprintf('lower_limit_ve = %s\n',num2str(prefs.lower_limit_ve));
    fprintf('upper_limit_ve = %s\n',num2str(prefs.upper_limit_ve));
    fprintf('initial_value_ve = %s\n',num2str(prefs.initial_value_ve));
    fprintf('lower_limit_tau = %s\n',num2str(prefs.lower_limit_tau));
    fprintf('upper_limit_tau = %s\n',num2str(prefs.upper_limit_tau));
    fprintf('initial_value_tau = %s\n',num2str(prefs.initial_value_tau));
    fprintf('TolFun = %s\n',num2str(prefs.TolFun));
    fprintf('TolX = %s\n',num2str(prefs.TolX));
    fprintf('MaxIter = %s\n',num2str(prefs.MaxIter));
    fprintf('MaxFunEvals = %s\n',num2str(prefs.MaxFunEvals));
    fprintf('Robust = %s\n',num2str(prefs.Robust));
    fprintf('fxr_fw = %s\n',num2str(prefs.fxr_fw));

    % Preallocate for speed
    GG = zeros([number_voxels 10],'double');
    % Slice out needed variables for speed
    Ct_data = xdata{1}.Ct;
    Cp_data = xdata{1}.Cp;
    timer_data = xdata{1}.timer;
    R1o= xdata{1}.R1o;
    R1i= xdata{1}.R1i;
    r1 = xdata{1}.relaxivity;
    fw = prefs.fxr_fw;
    %Turn off diary if on as it doesn't work with progress bar
    diary_restore = 0;
    if strcmp(get(0,'Diary'),'on')
        diary off;
        diary_restore = 1;
    end
    p = ProgressBar(number_voxels);
    parfor i = 1:number_voxels
        GG(i,:) = fxr_helper(Ct_data(:,i),Cp_data,timer_data,R1o(i),R1i(i),r1,fw,prefs);
        p.progress;
    end;
    p.stop;
    if diary_restore, diary on, end;
else
    warning(['Error, model ' model ' not yet implemented']);
    return
end


end
