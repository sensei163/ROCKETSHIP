function sizer = longestfilename(batchdata);
sizer = 0;

for i = 1:numel(batchdata)
    
    tempfiles = batchdata(i).files;
    
    for j = 1:numel(tempfiles)
        
        namer = tempfiles{j};
        [~, namer, ~] = fileparts(namer);
        sizer = max(sizer, numel(namer));
    end
end