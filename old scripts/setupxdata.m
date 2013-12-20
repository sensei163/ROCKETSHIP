%% b) ROI fitting
xdata{1}.timer = timer(starter:timeend)';
timer          = xdata{1}.timer;
CpROI          = mean(Cp,2);
CpROI          = CpROI(starter:timeend);
xdata{1}.Cp    = CpROI';

Cp1 = Cp1(starter:timeend,:);
%Cp2 = Cp2(starter:timeend,:);

%% threshold (Remove noise manually)
if(threshold)
    
    ind = find(CpROI > threshold);
    
    
    for j = 1:numel(ind)
        
        CpROI(ind(j)) = [];
        timer(ind(j)) = [];
        Cp1(ind(j), :)= [];
        Cp2(ind(j),:) = [];
    end
    
    xdata{1}.Cp    = CpROI';
    xdata{1}.timer = timer;
end

