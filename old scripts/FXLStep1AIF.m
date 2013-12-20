%% FXLStep1AAIF, no vp
function Ct = FXLStep1AIF(x,xdata)

Cp = xdata{1}.Cp;
T1 = xdata{1}.timer;

Cp = Cp(:);
T1 = T1(:);
% 
Ktrans = x(1);
ve     = x(2);

% Pre-alocate for speed
Ct = zeros(1,numel(T1));
for k = 1:numel(T1)
    
    % The time for T
    T = T1(1:k);
    CP= Cp(1:k);
    
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
    
