function cursubjectID = findsubjectID(filename);

for i = 1:numel(filename)
    
    curname = filename(i:end);
    
    testname = strrep(curname, '-', '_'); % Otherwise it becomes a negative sign
    testname = str2num(testname);
    
    if ~isempty(testname)
        cursubjectID = filename(1:(i-1));
        return;
    end
end