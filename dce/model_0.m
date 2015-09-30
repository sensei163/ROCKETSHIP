
function [x, residuals] = model_0(Ct)
  
% mean_ct = mean(Ct);

residuals = Ct - 0;
sse = sum(residuals.^2);

x(1) = sse;                 % residual

