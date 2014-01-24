function subjectID = finduniqueIDhelper(current, pre)

if isempty(pre)
    %initialize
    subjectID = current;
else
    
    %Search through the whole string of pre to find the longest contiguous
    %string that is found is current
    
    IDlength = 0;
    
    for i = 1:numel(pre)
        
        for j = i:numel(pre)
            
            compstr = pre(i:j);
            
            ans = strfind(current, compstr);
            
            if ~isempty(ans)
                ans = compstr;
            end
            
            if length(ans) > IDlength
                subjectID = ans;
                IDlength = length(ans);
            else
            end
        end
    end
end
            