function subsets = split2subsets(sort_list);

tempsort_list = sort_list;
subsettotal = 0;
subsets = [];

while size(tempsort_list,1) > 1
    
    cursubjectID = finduniqueIDhelper(tempsort_list{2}, tempsort_list{1});
    
    for i = 1:numel(tempsort_list)
        
        if ~isempty(strfind(tempsort_list{i}, cursubjectID))
            
            subsettotal = subsettotal +1;
        else
        end
    end
    
    subsets(end+1) = subsettotal;
    
    tempsort_list(1:subsettotal) = [];
    subsettotal = 0;
end
            
    