%% Smooth curves according to brix

function [sm x] = smoothcurve(curve, xdata, zeroed)
warning off
inj  = xdata{1}.inject;
timer=xdata{1}.timer;
oldtimer = timer;

ind = find(timer<= inj);
ind = ind(end);

%Set the preinjection value to 0 if concentration based, or R1
if(zeroed)
    xdata{1}.preinj = mean(curve(3:ind));
else
    xdata{1}.preinj = 0;
end

xdata{1}.preinj;

curve = curve(:);

ind = find(curve < mean(curve(3:ind)));
curve(ind) = [];
timer(ind) = [];
xdata{1}.timer = timer;
W = ones(size(curve));

if(numel(curve) < 0.4*numel(oldtimer))
    % Sparse, return nothing
    sm = -1;
    x = -1;
else
    
    ind = find(curve == max(curve(:)));
    
    W(ind) = 1;
    W(ind+1)= 1;
    
    
    scale = max(curve(:));
    
    curve = curve./scale;
    
    %configure the optimset for use with lsqcurvefit
    options = optimset('lsqcurvefit');
    
    %increase the number of function evaluations for more accuracy
    options.MaxFunEvals = 50;
    options.MaxIter     = 50;
    options.TolFun      = 10^(-20);
    options.TolX        = 10^(-20);
    options.Diagnostics = 'off';
    options.Display     = 'off';
    options.Algorithm   = 'levenberg-marquardt';
    
    [x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@brixfit, ...
        [1 0.5 0.9 inj], xdata, ...
        (curve)','','',options);
    rsquare = 1 - resnorm / norm(curve-mean(curve))^2;
    xdata{1}.timer = oldtimer;
    sm = brixfit(x, xdata)';
    
    sm = sm.*scale;
    
    ind = find(sm<xdata{1}.preinj);
    sm(ind) = xdata{1}.preinj;
    

    
   % figure, plot(timer, curve), hold on, plot(oldtimer, sm, 'g');
end


