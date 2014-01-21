function subjectID = finduniqueID(fullpath)

subjectID = '';
for i = 1:numel(fullpath)
    
    [~, filename, ~] = fileparts(fullpath{i});
    
    subjectID = finduniqueIDhelper(filename, subjectID);
end