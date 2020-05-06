function gt = gfun(t0,tmax,ymax,alpha,time)

%This function, GFUN, is used to define a fitting function. 

gt = zeros(1,numel(time)); 

for i = 1 : numel(time) 
    
    if time(i) < t0 
    gt(i) = ymax * (( 0 / (tmax - t0))^alpha) * exp( alpha - alpha * (0 / (tmax - t0))) ;   
    else
    gt(i) = ymax * (( (time(i) - t0) / (tmax - t0))^alpha) * exp( alpha - alpha * (time(i) - t0) / (tmax - t0)); 
    end
    
end
gt = gt';  
    
    
