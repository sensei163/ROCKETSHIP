%% FXLStep1AAIF, no vp
function Ct = FXLStep1AIFcfit(Ktrans, ve, Cp, T1)

% Cp = x{1}.Cp;
% T1 = x{1}.timer;

Cp = Cp(:);
T1 = T1(:);
% % 
%  Ktrans = a;
%  ve     = b;

% Pre-alocate for speed
Ct = zeros(1,numel(T1));
for k = 1:numel(T1)
    
    % The time for T
    T = T1(1:k);
    CP= Cp(1:k);
    
%     T = T(:);
%     CP= CP(:);
    
    F = CP.*exp((-Ktrans./ve).*(T(end)-T));
    
%     M = sampleintegration(T,F);
	if(numel(T) == 1)
		%need this as trapz interprets non array as
		%Y,DIM input instead of X,Y
		M = 0;
	else
		% 54 times faster than sampleintegration
		M = trapz(T,F);
	end
	
    Ct(k) = Ktrans.*M;
end
    
Ct = Ct';