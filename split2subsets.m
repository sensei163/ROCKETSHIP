function subsets = split2subsets(subjectID_current, rootname);

tempsort_list = subjectID_current;
subsettotal = 0;
subsets = [];



if numel(tempsort_list) == 1
    subsets(1) = 1;
% elseif allnumeric(tempsort_list)
%     subsets(1) = numel(tempsort_list);
%     % All numeric values, belong to the same volume
else
    
    while numel(tempsort_list) > 1
 
        pre = tempsort_list{1};
        post= tempsort_list{2};
        
        % Checks if the subject ID for each file is exactly the same as
        % rootname, if so, don't take out as there will be nothing to
        % compare.
        if ~strcmp(pre, rootname)
        pre = strrep(pre, rootname, '');
        end
        if ~strcmp(post, rootname)
        post= strrep(post, rootname, '');
        end
        
        cursubjectID = finduniqueIDhelper(post, pre);
        
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
end
            
    