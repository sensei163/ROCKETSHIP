function [GG, residuals] = FXLfit_generic(xdata, number_voxels, model, verbose)

if nargin < 4
    verbose = 1;
end
p = 0;
residuals = [];
if strcmp(model, 'ex_tofts')
    
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
    if verbose
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
    end
    
    % Preallocate for speed
    GG = zeros([number_voxels 10],'double');
    residuals = zeros([number_voxels numel(xdata{1}.timer)],'double');
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
    if verbose; p = ProgressBar(number_voxels,'verbose',verbose); end;
    parfor i = 1:number_voxels
        [GG(i,:), residuals(i,:)] = model_extended_tofts(Ct_data(:,i),Cp_data,timer_data,prefs);
        if verbose; p.progress; end;
    end;
    if verbose; p.stop; end;
    if diary_restore, diary on, end;
elseif strcmp(model, 'tissue_uptake')

    % Get values from pref file
    prefs_str = parse_preference_file('dce_preferences.txt',0,...
        {'voxel_lower_limit_ktrans' 'voxel_upper_limit_ktrans' 'voxel_initial_value_ktrans' ...
        'voxel_lower_limit_fp' 'voxel_upper_limit_fp' 'voxel_initial_value_fp' ...
        'voxel_lower_limit_tp' 'voxel_upper_limit_tp' 'voxel_initial_value_tp' ...
        'voxel_TolFun' 'voxel_TolX' 'voxel_MaxIter' 'voxel_MaxFunEvals' 'voxel_Robust'});
    prefs.lower_limit_ktrans = str2num(prefs_str.voxel_lower_limit_ktrans);
    prefs.upper_limit_ktrans = str2num(prefs_str.voxel_upper_limit_ktrans);
    prefs.initial_value_ktrans = str2num(prefs_str.voxel_initial_value_ktrans);
    prefs.lower_limit_fp = str2num(prefs_str.voxel_lower_limit_fp);
    prefs.upper_limit_fp = str2num(prefs_str.voxel_upper_limit_fp);
    prefs.initial_value_fp = str2num(prefs_str.voxel_initial_value_fp);
    prefs.lower_limit_tp = str2num(prefs_str.voxel_lower_limit_tp);
    prefs.upper_limit_tp = str2num(prefs_str.voxel_upper_limit_tp);
    prefs.initial_value_tp = str2num(prefs_str.voxel_initial_value_tp);
    prefs.TolFun = str2num(prefs_str.voxel_TolFun);
    prefs.TolX = str2num(prefs_str.voxel_TolX);
    prefs.MaxIter = str2num(prefs_str.voxel_MaxIter);
    prefs.MaxFunEvals = str2num(prefs_str.voxel_MaxFunEvals);
    prefs.Robust = prefs_str.voxel_Robust;
    %Log values used
    if verbose
        fprintf('lower_limit_ktrans = %s\n',num2str(prefs.lower_limit_ktrans));
        fprintf('upper_limit_ktrans = %s\n',num2str(prefs.upper_limit_ktrans));
        fprintf('initial_value_ktrans = %s\n',num2str(prefs.initial_value_ktrans));
        fprintf('lower_limit_fp = %s\n',num2str(prefs.lower_limit_fp));
        fprintf('upper_limit_fp = %s\n',num2str(prefs.upper_limit_fp));
        fprintf('initial_value_fp = %s\n',num2str(prefs.initial_value_fp));
        fprintf('lower_limit_tp = %s\n',num2str(prefs.lower_limit_tp));
        fprintf('upper_limit_tp = %s\n',num2str(prefs.upper_limit_tp));
        fprintf('initial_value_tp = %s\n',num2str(prefs.initial_value_tp));
        fprintf('TolFun = %s\n',num2str(prefs.TolFun));
        fprintf('TolX = %s\n',num2str(prefs.TolX));
        fprintf('MaxIter = %s\n',num2str(prefs.MaxIter));
        fprintf('MaxFunEvals = %s\n',num2str(prefs.MaxFunEvals));
        fprintf('Robust = %s\n',num2str(prefs.Robust));
    end
    
    % Preallocate for speed
    GG = zeros([number_voxels 10],'double');
    residuals = zeros([number_voxels numel(xdata{1}.timer)],'double');
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
    if verbose; p = ProgressBar(number_voxels,'verbose',verbose); end;
    parfor i = 1:number_voxels
        % Do quick linear patlak and use values as initial values
        [estimate, ~] = model_patlak_linear(Ct_data(:,i),Cp_data,timer_data);
        prefs_local = prefs;
        prefs_local.initial_value_ktrans = estimate(1);
        prefs_local.initial_value_vp = estimate(2);
        % Do tissue uptake fit
        [GG(i,:), residuals(i,:)] = model_tissue_uptake(Ct_data(:,i),Cp_data,timer_data,prefs_local);
        if verbose; p.progress; end;
    end;
    if verbose; p.stop; end;
    if diary_restore, diary on, end;
elseif strcmp(model, 'tofts')
    
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
    if verbose
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
    end
    
    % Preallocate for speed
    GG = zeros([number_voxels 7],'double');
    residuals = zeros([number_voxels numel(xdata{1}.timer)],'double');
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
    if verbose; p = ProgressBar(number_voxels,'verbose',verbose); end;
    parfor i = 1:number_voxels
        [GG(i,:), residuals(i,:)] = model_tofts(Ct_data(:,i),Cp_data,timer_data,prefs);
        if verbose; p.progress; end;
    end;
    if verbose; p.stop; end;
    if diary_restore, diary on, end;
elseif strcmp(model, 'fxr')
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
    if verbose
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
    end
    
    % Preallocate for speed
    GG = zeros([number_voxels 10],'double');
    residuals = zeros([number_voxels numel(xdata{1}.timer)],'double');
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
    if verbose; p = ProgressBar(number_voxels,'verbose',verbose); end;
    parfor i = 1:number_voxels
        [GG(i,:), residuals(i,:)] = model_fxr(Ct_data(:,i),Cp_data,timer_data,R1o(i),R1i(i),r1,fw,prefs);
        if verbose; p.progress; end;
    end;
    if verbose; p.stop; end;
    if diary_restore, diary on, end;
elseif strcmp(model, 'auc')
    % Area under curve using Raw data signal
    
    % Preallocate for speed
    GG = zeros([number_voxels 4],'double');
    % Slice out needed variables for speed
    Sss    = xdata{1}.Sss;
    Ssstum = xdata{1}.Ssstum;
    Stlv   = xdata{1}.Stlv;
    Sttum  = xdata{1}.Sttum;
    Ct_data = xdata{1}.Ct;
    Cp_data = xdata{1}.Cp;
    
    % Substract steady state signal from time curves
    Sttum = Sttum - repmat(Ssstum, [size(Sttum,1) 1]);
    Sss = mean(Sss);
    Stlv  = Stlv - Sss;%repmat(Sss, [size(Stlv,1) 1]);
    %Stlv  = mean(Stlv, 2);
    
    timer_data      = xdata{1}.timer;
    start_injection = xdata{1}.start_injection;
    end_injection   = xdata{1}.end_injection;
    
    % Extract signal only after injection has started
    ind = find(timer_data >= start_injection);
    
    Sttum      = Sttum(ind(1):end,:);
    Stlv       = Stlv(ind(1):end);
    Ct_data    = Ct_data(ind(1):end,:);
    Cp_data    = Cp_data(ind(1):end);
    timer_data = timer_data(ind(1):end);
    
    %Turn off diary if on as it doesn't work with progress bar
    diary_restore = 0;
    if strcmp(get(0,'Diary'),'on')
        diary off;
        diary_restore = 1;
    end
    if verbose; p = ProgressBar(number_voxels,'verbose',verbose); end;
    parfor i = 1:number_voxels
        GG(i,:) = auc_helper(Sttum(:,i),Stlv, Ct_data(:,i), Cp_data, timer_data);
        if verbose; p.progress; end;
    end;
    if verbose; p.stop; end;
    if diary_restore, diary on, end;
    
elseif strcmp(model, 'fractal')
    
elseif strcmp(model, 'nested')
    
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
    if verbose
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
    end
    
    % Preallocate for speed
    GG = zeros([number_voxels 10],'double');
    residuals = zeros([number_voxels numel(xdata{1}.timer)],'double');
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
    if verbose; p = ProgressBar(number_voxels,'verbose',verbose); end;
    
    parfor i = 1:number_voxels
        % Fit 0 order model
        [GG_zero, residuals(i,:)] = model_0(Ct_data(:,i));
        fp_lower = 0;
        % Fit 1st order model
        [GG_one, residuals_b] = model_vp(Ct_data(:,i),Cp_data,timer_data,prefs);
        fp_higher = 1;
        % Compare
        n = numel(timer_data);
        p_value = ftest(n,fp_lower,fp_higher,GG_zero,GG_one(1,2));
        if p_value>=0.05
            % Use 0 order
            GG_local = zeros([1 10],'double');
            GG_local(1,4) = GG_zero;
            GG(i,:) = GG_local;
            residuals(i,:) = residuals_b';
        else
            % Use 1st order
            GG_local = zeros([1 10],'double');
            GG_local(1,3) = GG_one(1,1);
            GG_local(1,4) = GG_one(1,2);
            GG_local(1,9) = GG_one(1,3);
            GG_local(1,10) = GG_one(1,4);
            GG(i,:) = GG_local;
            residuals(i,:) = residuals_b';
            fp_lower = 1;
            
            % Continue, Fit 2nd order model
            [GG_two, residuals_b] = model_patlak(Ct_data(:,i),Cp_data,timer_data,prefs);
            fp_higher = 2;
            % Compare
            p_value = ftest(n,fp_lower,fp_higher,GG_one(1,2),GG_two(1,3));
            if p_value<0.05
                % Use 2nd order
                GG_local = zeros([1 10],'double');
                GG_local(1,1) = GG_two(1,1);
                GG_local(1,3) = GG_two(1,2);
                GG_local(1,4) = GG_two(1,3);
                GG_local(1,5) = GG_two(1,4);
                GG_local(1,6) = GG_two(1,5);
                GG_local(1,9) = GG_two(1,6);
                GG_local(1,10) = GG_two(1,7);
                GG(i,:) = GG_local;
                residuals(i,:) = residuals_b';
                fp_lower = 2;
                
                % Continue, Fit 3rd order model
                [GG_three, residuals_b] = model_extended_tofts(Ct_data(:,i),Cp_data,timer_data,prefs);
                fp_higher = 3;
                % Compare
                p_value = ftest(n,fp_lower,fp_higher,GG_two(1,3),GG_three(1,4));
                if p_value<0.05
                    % Use 3rd order
                    GG(i,:) = GG_three;
                    residuals(i,:) = residuals_b';
                    fp_lower = 3;
                end
            end
        end
        
        if verbose; p.progress; end;
    end;
    if verbose; p.stop; end;
    if diary_restore, diary on, end;
    
elseif strcmp(model, 'patlak')
    
    % Get values from pref file
    prefs_str = parse_preference_file('dce_preferences.txt',0,...
        {'voxel_lower_limit_ktrans' 'voxel_upper_limit_ktrans' 'voxel_initial_value_ktrans' ...
        'voxel_lower_limit_vp' 'voxel_upper_limit_vp' 'voxel_initial_value_vp' ...
        'voxel_TolFun' 'voxel_TolX' 'voxel_MaxIter' 'voxel_MaxFunEvals' 'voxel_Robust'});
    prefs.lower_limit_ktrans = str2num(prefs_str.voxel_lower_limit_ktrans);
    prefs.upper_limit_ktrans = str2num(prefs_str.voxel_upper_limit_ktrans);
    prefs.initial_value_ktrans = str2num(prefs_str.voxel_initial_value_ktrans);
    prefs.lower_limit_vp = str2num(prefs_str.voxel_lower_limit_vp);
    prefs.upper_limit_vp = str2num(prefs_str.voxel_upper_limit_vp);
    prefs.initial_value_vp = str2num(prefs_str.voxel_initial_value_vp);
    prefs.TolFun = str2num(prefs_str.voxel_TolFun);
    prefs.TolX = str2num(prefs_str.voxel_TolX);
    prefs.MaxIter = str2num(prefs_str.voxel_MaxIter);
    prefs.MaxFunEvals = str2num(prefs_str.voxel_MaxFunEvals);
    prefs.Robust = prefs_str.voxel_Robust;
    %Log values used
    if verbose
        fprintf('lower_limit_ktrans = %s\n',num2str(prefs.lower_limit_ktrans));
        fprintf('upper_limit_ktrans = %s\n',num2str(prefs.upper_limit_ktrans));
        fprintf('initial_value_ktrans = %s\n',num2str(prefs.initial_value_ktrans));
        fprintf('lower_limit_vp = %s\n',num2str(prefs.lower_limit_vp));
        fprintf('upper_limit_vp = %s\n',num2str(prefs.upper_limit_vp));
        fprintf('initial_value_vp = %s\n',num2str(prefs.initial_value_vp));
        fprintf('TolFun = %s\n',num2str(prefs.TolFun));
        fprintf('TolX = %s\n',num2str(prefs.TolX));
        fprintf('MaxIter = %s\n',num2str(prefs.MaxIter));
        fprintf('MaxFunEvals = %s\n',num2str(prefs.MaxFunEvals));
        fprintf('Robust = %s\n',num2str(prefs.Robust));
    end
    
    % Preallocate for speed
    GG = zeros([number_voxels 7],'double');
    residuals = zeros([number_voxels numel(xdata{1}.timer)],'double');
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
    if verbose; p = ProgressBar(number_voxels,'verbose',verbose); end;
    parfor i = 1:number_voxels
        % Do quick linear patlak and use values as initial values
        [estimate, ~] = model_patlak_linear(Ct_data(:,i),Cp_data,timer_data);
        prefs_local = prefs;
        prefs_local.initial_value_ktrans = estimate(1);
        prefs_local.initial_value_vp = estimate(2);
        % Do non-linear patlak
        [GG(i,:), residuals(i,:)] = model_patlak(Ct_data(:,i),Cp_data,timer_data,prefs_local);
        if verbose; p.progress; end;
    end;
    if verbose; p.stop; end;
    if diary_restore, diary on, end;

elseif strcmp(model, 'patlak_linear')

    % Preallocate for speed
    GG = zeros([number_voxels 7],'double');
    residuals = zeros([number_voxels numel(xdata{1}.timer)],'double');
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
    if verbose; p = ProgressBar(number_voxels,'verbose',verbose); end;
    parfor i = 1:number_voxels
        [GG(i,:), residuals(i,:)] = model_patlak_linear(Ct_data(:,i),Cp_data,timer_data);
        if verbose; p.progress; end;
    end;
    if verbose; p.stop; end;
    if diary_restore, diary on, end;
    
elseif strcmp(model, '2cxm')
    
    % Get values from pref file
    prefs_str = parse_preference_file('dce_preferences.txt',0,...
        {'voxel_lower_limit_ktrans' 'voxel_upper_limit_ktrans' 'voxel_initial_value_ktrans' ...
        'voxel_lower_limit_ve' 'voxel_upper_limit_ve' 'voxel_initial_value_ve' ...
        'voxel_lower_limit_vp' 'voxel_upper_limit_vp' 'voxel_initial_value_vp' ...
        'voxel_lower_limit_fp' 'voxel_upper_limit_fp' 'voxel_initial_value_fp' ...
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
    prefs.lower_limit_fp = str2num(prefs_str.voxel_lower_limit_fp);
    prefs.upper_limit_fp = str2num(prefs_str.voxel_upper_limit_fp);
    prefs.initial_value_fp = str2num(prefs_str.voxel_initial_value_fp);
    prefs.TolFun = str2num(prefs_str.voxel_TolFun);
    prefs.TolX = str2num(prefs_str.voxel_TolX);
    prefs.MaxIter = str2num(prefs_str.voxel_MaxIter);
    prefs.MaxFunEvals = str2num(prefs_str.voxel_MaxFunEvals);
    prefs.Robust = prefs_str.voxel_Robust;
    %Log values used
    if verbose
        fprintf('lower_limit_ktrans = %s\n',num2str(prefs.lower_limit_ktrans));
        fprintf('upper_limit_ktrans = %s\n',num2str(prefs.upper_limit_ktrans));
        fprintf('initial_value_ktrans = %s\n',num2str(prefs.initial_value_ktrans));
        fprintf('lower_limit_ve = %s\n',num2str(prefs.lower_limit_ve));
        fprintf('upper_limit_ve = %s\n',num2str(prefs.upper_limit_ve));
        fprintf('initial_value_ve = %s\n',num2str(prefs.initial_value_ve));
        fprintf('lower_limit_vp = %s\n',num2str(prefs.lower_limit_vp));
        fprintf('upper_limit_vp = %s\n',num2str(prefs.upper_limit_vp));
        fprintf('initial_value_vp = %s\n',num2str(prefs.initial_value_vp));
        fprintf('lower_limit_fp = %s\n',num2str(prefs.lower_limit_fp));
        fprintf('upper_limit_fp = %s\n',num2str(prefs.upper_limit_fp));
        fprintf('initial_value_fp = %s\n',num2str(prefs.initial_value_fp));
        fprintf('TolFun = %s\n',num2str(prefs.TolFun));
        fprintf('TolX = %s\n',num2str(prefs.TolX));
        fprintf('MaxIter = %s\n',num2str(prefs.MaxIter));
        fprintf('MaxFunEvals = %s\n',num2str(prefs.MaxFunEvals));
        fprintf('Robust = %s\n',num2str(prefs.Robust));
    end
    
    % Preallocate for speed
    GG = zeros([number_voxels 13],'double');
    residuals = zeros([number_voxels numel(xdata{1}.timer)],'double');
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
    if verbose; p = ProgressBar(number_voxels,'verbose',verbose); end;
    parfor i = 1:number_voxels
%         [estimate, ~] = model_patlak_linear(Ct_data(:,i),Cp_data,timer_data);
%         prefs_local = prefs;
%         prefs_local.initial_value_ktrans = estimate(1);
        [GG(i,:), residuals(i,:)] = model_2cxm(Ct_data(:,i),Cp_data,timer_data,prefs);
        if verbose; p.progress; end;
    end;
    if verbose; p.stop; end;
    if diary_restore, diary on, end;
    
elseif strcmp(model, 'FXL_rr')
    
    % Get values from pref file
    prefs_str = parse_preference_file('dce_preferences.txt',0,...
        {'voxel_lower_limit_ktrans' 'voxel_upper_limit_ktrans' 'voxel_initial_value_ktrans' ...
        'voxel_lower_limit_ve' 'voxel_upper_limit_ve' 'voxel_initial_value_ve' ...
        'voxel_lower_limit_ktrans_RR' 'voxel_upper_limit_ktrans_RR' 'voxel_initial_value_ktrans_RR' ...
        'voxel_value_ve_RR' ...
        'voxel_TolFun' 'voxel_TolX' 'voxel_MaxIter' 'voxel_MaxFunEvals' 'voxel_Robust'});
    prefs.lower_limit_ktrans = str2num(prefs_str.voxel_lower_limit_ktrans);
    prefs.upper_limit_ktrans = str2num(prefs_str.voxel_upper_limit_ktrans);
    prefs.initial_value_ktrans = str2num(prefs_str.voxel_initial_value_ktrans);
    prefs.lower_limit_ve = str2num(prefs_str.voxel_lower_limit_ve);
    prefs.upper_limit_ve = str2num(prefs_str.voxel_upper_limit_ve);
    prefs.initial_value_ve = str2num(prefs_str.voxel_initial_value_ve);
    prefs.lower_limit_ktrans_RR = str2num(prefs_str.voxel_lower_limit_ktrans_RR);
    prefs.upper_limit_ktrans_RR = str2num(prefs_str.voxel_upper_limit_ktrans_RR);
    prefs.initial_value_ktrans_RR = str2num(prefs_str.voxel_initial_value_ktrans_RR);
    prefs.value_ve_RR = str2num(prefs_str.voxel_value_ve_RR);
    prefs.TolFun = str2num(prefs_str.voxel_TolFun);
    prefs.TolX = str2num(prefs_str.voxel_TolX);
    prefs.MaxIter = str2num(prefs_str.voxel_MaxIter);
    prefs.MaxFunEvals = str2num(prefs_str.voxel_MaxFunEvals);
    prefs.Robust = prefs_str.voxel_Robust;
    %Log values used
    if verbose
        fprintf('lower_limit_ktrans = %s\n',num2str(prefs.lower_limit_ktrans));
        fprintf('upper_limit_ktrans = %s\n',num2str(prefs.upper_limit_ktrans));
        fprintf('initial_value_ktrans = %s\n',num2str(prefs.initial_value_ktrans));
        fprintf('lower_limit_ve = %s\n',num2str(prefs.lower_limit_ve));
        fprintf('upper_limit_ve = %s\n',num2str(prefs.upper_limit_ve));
        fprintf('initial_value_ve = %s\n',num2str(prefs.initial_value_ve));
        fprintf('lower_limit_ktrans_RR = %s\n',num2str(prefs.lower_limit_ktrans_RR));
        fprintf('upper_limit_ktrans_RR = %s\n',num2str(prefs.upper_limit_ktrans_RR));
        fprintf('initial_value_ktrans_RR = %s\n',num2str(prefs.initial_value_ktrans_RR));
        fprintf('value_ve_RR = %s\n',num2str(prefs.value_ve_RR));
        fprintf('TolFun = %s\n',num2str(prefs.TolFun));
        fprintf('TolX = %s\n',num2str(prefs.TolX));
        fprintf('MaxIter = %s\n',num2str(prefs.MaxIter));
        fprintf('MaxFunEvals = %s\n',num2str(prefs.MaxFunEvals));
        fprintf('Robust = %s\n',num2str(prefs.Robust));
    end
    
    % Preallocate for speed
    GG = zeros([number_voxels 10],'double');
    residuals = zeros([number_voxels numel(xdata{1}.timer)],'double');
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
    if verbose; p = ProgressBar(number_voxels,'verbose',verbose); end;
    parfor i = 1:number_voxels 

        [GG(i,:), residuals(i,:)] = model_FXL_reference_region(Ct_data(:,i),Cp_data,timer_data,prefs);
        if verbose; p.progress; end;
    end;
    if verbose; p.stop; end;
    if diary_restore, diary on, end;
    
else
    error(['Error, model ' model ' not yet implemented']);
    return
end
