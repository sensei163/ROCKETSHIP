% Remove dir info to allow nice visualization, DCE version
function [visual_listA, visual_listB, visual_listC] = visualize_list_dce(handles, filevolume)

visual_listA = [];
visual_listB = [];
visual_listC = [];

batchdata = handles.batchdata;

if filevolume == 1
    namer = batchdata(1).files;
    [~, namer, ext] = fileparts(namer{1});
    visual_listA = namer;
    
elseif filevolume == 2
    
    sizer = longestfilename(batchdata);
    
    for i = 1:numel(batchdata)
        
        namer = batchdata(i).files;
        [~, namer, ext] = fileparts(namer{1});
        
        if numel(namer) < sizer
            namer = [namer blanks(sizer-numel(namer))];
        end
        visual_listB = [visual_listB; namer];
        visual_listA = [visual_listA; namer];
    end
elseif filevolume == 3
    sizer = longestfilename(batchdata);
    
    for i = 1:numel(batchdata)
        
        files = batchdata{i}.files;
        for j = 1:numel(files)
            namer = files{j};
            [~, namer, ext] = fileparts(namer);
            
            if numel(namer) < sizer
                namer = [namer blanks(sizer-numel(namer))];
            end
            
            if j == 1
                maxpad = numel(num2str(numel(batchdata)*numel(files)));
                pad = [blanks(maxpad-numel(num2str(numel(files)))) num2str(numel(files))];
                namerB = [namer ' --> ' pad ' files'];
                visual_listB = [visual_listB; namerB];
            end
            
            visual_listC = [visual_listC; namer];
        end
    end
else
    disp('invalid file volume defined');
end


