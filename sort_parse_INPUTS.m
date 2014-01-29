function [sortedfullpath, sortlist_current, subjectID_current] = sort_parse_INPUTS(fullpath, natdigits);

for i = 1:numel(fullpath)
    curpath = fullpath{i};
    
    [~, curpath, ~] = fileparts(curpath);
    
    cursubjectID = findsubjectID(curpath);
    
    oldsubjectID_current{i} = cursubjectID;
    curpath = strrep(curpath, cursubjectID, '');
    
    %Find all the numbers in the curpath
    
    numbers = regexp(curpath,['\d+'], 'match');
    
    placeholder = num2str(round(rand(1)*10^natdigits));
    
    tempstr = regexprep(curpath,['\d+'], placeholder);
    
    ind     = strfind(tempstr, placeholder);
    
    % Replace numbers
    for j = 1:numel(numbers)
        
        newnum = [num2str(zeros(1, natdigits-numel(numbers{j}))) numbers{j}];
        newnum = strrep(newnum,' ', '');
        
        tempstr(ind(j):(ind(j)+natdigits-1)) = newnum;
    end
    
    editlist{i} = [cursubjectID tempstr];
end

% Now sort the list

[sortlist_current, ind] = sort(editlist);%, 'ascend');

for i = 1:numel(fullpath)
    sortedfullpath(i)   = fullpath(ind(i));
    subjectID_current(i)= oldsubjectID_current(ind(i));
end

