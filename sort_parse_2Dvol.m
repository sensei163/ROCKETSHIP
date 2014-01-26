function [subsets, errormsgA] = sort_parse_2Dvol(subjectID_current, rootname);

errormsgA = '';

subsets = split2subsets(subjectID_current,  rootname);

% Check if there are uneven subsets

if numel(unique(subsets)) > 1
    errormsg = 'Uneven subsets: Check filenames';
    return;
end