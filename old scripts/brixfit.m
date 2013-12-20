%% Implements the Brix et al signal enhancement equation

function S = brixfit(x, xdata)

C1 = x(1);
C2 = x(2);
C3 = x(3);
C4 = x(4);

inject = xdata{1}.inject;
timer  = xdata{1}.timer;
preinj = xdata{1}.preinj;



for i = 1:numel(timer)
    t = timer(i);
    
    if(t <= inject)
        S(i) = preinj;
    else
        S(i) = (C1/(C2-C3))*(exp(-C3*(t-C4)) - exp(-C2*(t-C4)));
    end
    
end

