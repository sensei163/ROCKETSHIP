function cursubjectID = findsubjectID(filename);

% cursubjectID = findsubjectID(filename)
% findsubjectID attempts to separate the filename between the image counter
% and the "root" of the file. Outputs the root.

% We remove any lettering at the end of the filename first, assume that the
% files are ordered by numbers.

while(isempty(str2num(filename(end))))
    
    ind = strfind(filename, filename(end));
    filename(ind(end)) = '';
end
% No we sort through the rest of the file name

for i = 1:numel(filename)
    
    curname = filename(i:end);
    testname = strrep(curname, '-', '_'); % Otherwise it becomes a negative sign
    testname = str2num(testname);
    
    if ~isempty(testname)
        cursubjectID = filename(1:(i-1));
        return;
    end
end