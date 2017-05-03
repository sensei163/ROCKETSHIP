%% FXL Reference Region model

function Ct = model_FXL_reference_region_cfit(Ktrans_TOI, ve_TOI, Ktrans_RR, Cp, T1)

Cp   = Cp(:);
T1   = T1(:);
file_struct = parse_preference_file('dce_preferences.txt',0,...
    { 'voxel_value_ve_RR'});
ve_RR = str2num(file_struct.voxel_value_ve_RR);

%Pre allocate for speed

Ct = zeros(1,numel(T1));

R = Ktrans_TOI/Ktrans_RR;

for k = 1:numel(T1)
    
        T = T1(1:k);
        CP = Cp(1:k);
        
        F = CP.*(exp(-Ktrans_TOI./ve_TOI).*(T(end)-T));
               
        if(numel(T) == 1)
            M = 0;
        else
            M = trapz(T,F);
        end
        
        Ct(k) = R.*Cp(k) + R.*((Ktrans_RR./ve_RR) - Ktrans_TOI./ve_TOI).*M;
  
end

Ct = Ct';