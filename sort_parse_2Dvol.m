function [subjectIDA, handles, errormsg] = sort_parse_2Dvol(fullpath, handles, rootnameA);

if isempty(rootnameA)
    subjectIDA = finduniqueID(fullpath);
else
    subjectIDA = rootnameA;
end

% Sort all files first, just like with 3D volume

%number of digits to have in order to facilitate natural
%order sorting
natdigits = 6;

[sortedfullpath, sort_list] = sort_parse_3Dvol_helper(subjectID, fullpath, natdigits);

[sortedfullpath, errormsg] = sort_parse_2Dvol_helper(sortedfullpath, sort_list);

%Either time or slices. We sort by first number after the subject ID

 handles.batchdata = sortedfullpath;

 


