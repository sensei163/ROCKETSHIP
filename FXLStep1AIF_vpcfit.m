%% FXLStep1AAIF, vp
function Ct = FXLStep1AIF_vpcfit(Ktrans, ve, vp, Cp, T1)

% Cp = xdata{1}.Cp;
% T1 = xdata{1}.timer;

Cp = Cp(:);
T1 = T1(:);

% Ktrans = x(1);
% ve     = x(2);
% vp     = x(3);


for k = 1:numel(T1)
    
    % The time for T
    
    T = T1(1:k);
    
    CP= Cp(1:k);
    
    T = T(:);
    CP= CP(:);
    
    F = CP.*exp((-Ktrans./ve).*(T(end)-T));
    
    M = sampleintegration(T,F);
    
    Ct(k) = Ktrans.*M+vp*Cp(k);
end

Ct = Ct';
    