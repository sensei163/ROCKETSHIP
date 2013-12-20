%% b) Time management
% Convert time to minutes
timereso = timeres;
timeres = timeres/60;

%Generate time curve
timer   = [0:timeres:timeres*(size(Cp,1)-1)];

% %Starting time set (in minutes, 0 if all
% starter = 0;
% 
% %Time limit (in minutes), 0 if all
% timeend = 0;

if(starter)
    
    timecheck = abs(timer-starter);
    starter = find(timecheck == min(timecheck));
else
   % timeend = numel(timer);
    starter = 1;
end

%
if(timeend)
    
    timecheck = abs(timer-timeend);
    timeend = find(timecheck == min(timecheck));
else
    timeend = numel(timer);
end