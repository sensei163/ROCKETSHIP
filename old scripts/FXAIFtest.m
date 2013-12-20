%% FXR AIF model

function R1t = FXRAIFtest(Ktrans, ve, tau1, xdata,Cp, T1)
R1o = xdata{1}.R10;
R1i = xdata{1}.R1i;

fw  = xdata{1}.fw;
r1  = xdata{1}.r1;


po(1) = ve(1)/fw(1);

% Calculate Ct first
Cp = Cp(:);
T1 = T1(:);

for k = 1:numel(T1)
    
    % The time for T
    
    T = T1(1:k);
    
    CP= Cp(1:k);
    
    T = T(:);
    CP= CP(:);
    
    F = CP.*exp((-Ktrans./ve).*(T(end)-T));
    
    M = sampleintegration(T,F);
    
    Ct(k) = Ktrans.*M;
end
    
Ct = Ct';


%% Now incorporate into FXR 2SX Equation

XX = (R1o-R1i+1/tau1)/po;
YY = ((2/tau1 -r1.*Ct-XX).^2)+4*(1-po)/(tau1*tau1*po);

ind = find(YY < 0);

if(~isempty(ind))
    YY(ind) = 0;
else
end
R1t = (1/2).*(2.*R1i+r1.*Ct+XX-sqrt(YY));
