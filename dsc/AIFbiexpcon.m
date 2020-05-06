function Cp = AIFbiexpcon( x, xdata )
%% AIFbiexpconvolve Cp

%% This models the AIF as a biexponential function convolved with a
%% rectangular function to model the injection
A = x(1);
B = x(2);
c = x(3);
d = x(4);

baseline = xdata.baseline;
T1 = xdata.timer;
step = xdata.step;
bolus_time = xdata.bolus_time;
Cp_old = xdata.Cp;

OUT = find(step > 0);
max_index = OUT(end);



%Set as model A, MacGrath, MRM 2009

%B = maxer-A;

%{
if isfield(xdata{1}, 'raw') && xdata{1}.raw == true && isfield(xdata{1}, 'baseline')
    baseline = xdata{1}.baseline;
else
    baseline = 0;
end
%}

%find the index of the maxiumum value to control linear increase of the
%fitting function 
%[max_value, max_index] = max(Cp_old);
for j = 1:numel(T1)
    % Baseline
    if(j< bolus_time)
        Cp(j) = baseline;
    % Linear upslope to max
    elseif(j <= (max_index))
        Cp(j) = baseline + (A-baseline).*(j-bolus_time)./(max_index-bolus_time) + (B-baseline).*(j-bolus_time)./(max_index-bolus_time);
    % Bi-Exponential Decay    
    else
        Cp(j) = A.*(exp(-c.*(j - max_index))) + B.*(exp(-d.*(j - max_index)));
    end
end
end

