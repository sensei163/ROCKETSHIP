function [TUMOR, LV, NOISE, T1MAP, DYNAMIC, hdr, res, errormsg] = loadIMGVOL(handles);

% Takes handles, loads the image files and outputs image volume.

filevolume = get(handles.filevolume, 'Value');

noise_pathpick= get(handles.noisefile, 'Value');
noise_path = get(handles.noise_path,'String');
noise_pixelspick=get(handles.noisepixels, 'Value');
noise_pixsize=str2num(get(handles.noisepixsize, 'String'));

LUT = handles.LUT;

t1aiffiles = handles.t1aiffiles;
t1roifiles = handles.t1roifiles;
t1mapfiles = handles.t1mapfiles;
noisefiles = handles.noisefiles;
filelist   = handles.filelist;

fileorder = get(get(handles.fileorder,'SelectedObject'),'Tag');

quant     = get(handles.quant, 'Value');

% Initialize Image sets
TUMOR = [];
LV    = [];
NOISE = [];
T1MAP = [];
DYNAMIC=[];

%% Load image files

% Load AIF - either 3D volume or 2D slice
for i = 1:numel(t1aiffiles)
    
    if isDICOM(t1aiffiles{i})
        hdr = dicominfo(t1aiffiles{i});
        img = dicomread(hdr);
        LV  = img;
        
        
    elseif isNIFTI(t1aiffiles{i})
        nii = load_untouch_nii(t1aiffiles{i});
        img = nii.img;
        if i == 1
            LV = img;
        else
            LV(:,:,end+1) = img;
        end
    else
        errormsg = 'Unknown file type - AIF';
        return;
    end
end

% Load ROI file - either 3D volume or 2D slice
% hdr , res are derived from here

for i = 1:numel(t1roifiles)
    
    if isDICOM(t1roifiles{i})
        hdr = dicominfo(t1roifiles{i});
        img = dicomread(hdr);
        res(1:2) = hdr.(dicomlookup('28', '30'));
        res(3)   = hdr.(dicomlookup('18', '50'));
        TUMOR= img;
        
        
    elseif isNIFTI(t1roifiles{i})
        nii = load_untouch_nii(t1roifiles{i});
        hdr = nii.hdr;
        res = nii.hdr.dime.pixdim;
        res = res(2:4);
        img = nii.img;
        if i == 1
            TUMOR = img;
        else
            TUMOR(:,:,end+1) = img;
        end
    else
        errormsg = 'Unknown file type - ROI';
        return;
    end
end

if quant
    % Load T1 map if quantitative parameters desired
    for i = 1:numel(t1mapfiles)
        if isDICOM(t1mapfiles{i})
            hdr = dicominfo(t1mapfiles{i});
            img = dicomread(hdr);
            T1MAP  = img;
            
            
        elseif isNIFTI(t1aiffiles{i})
            nii = load_untouch_nii(t1aiffiles{i});
            img = nii.img;
            if i == 1
                T1MAP = img;
            else
                T1MAP(:,:,end+1) = img;
            end
        else
            errormsg = 'Unknown file type - AIF';
            return;
        end
    end
end

if noise_pathpick
    % Load noise map
    
    for i = 1:numel(noisefiles)
        if isDICOM(noisefiles{i})
            hdr = dicominfo(noisefiles{i});
            img = dicomread(hdr);
            NOISE  = img;
            
            
        elseif isNIFTI(noisefiles{i})
            nii = load_untouch_nii(noisefiles{i});
            img = nii.img;
            if i == 1
                NOISE = img;
            else
                NOISE(:,:,end+1) = img;
            end
        else
            errormsg = 'Unknown file type - NOISE';
            return;
        end
    end
end

%% Check to make sure the dimensions of the files are consistent

if ~isequal(size(TUMOR), size(LV)) || ~isequal(size(TUMOR), size(NOISE)) || ~isequal(size(TUMOR), size(T1MAP)) ...
        || ~isequal(size(LV), size(NOISE)) || ~isequal(size(LV), size(T1MAP)) || ~isequal(size(NOISE, T1MAP))
    errormsg = 'Input medical sizes not the same';
    return;
end
%% Now if the noise file is not defined, we make one from the pixelsize

if ~noise_pathpick
    
    NOISE = zeros(size(TUMOR));
    
    for i = 1:size(NOISE,3)
        
        NOISE(1:noise_pixsize, 1:noise_pixsize, i) = 1;
    end
end

%% Now we load dynamic

if filevolume == 1
    %4D
    id = LUT(1,1);
    
    if isDICOM(filelist{id})
        hdr = dicominfo(filelist{id});
        img = dicomread(hdr);
        DYNAMIC  = img;
        
        
    elseif isNIFTI(filelist{id})
        nii = load_untouch_nii(filelist{id});
        img = nii.img;
        
        DYNAMIC = img;
    else
        errormsg = 'Unknown file type - DYNAMIC';
        return;
    end
elseif filevolume == 2
    %3D
    
    if size(LUT,1) == 1
        for i = 1:size(LUT,2)
            id = LUT(1,i);
            if isDICOM(filelist{id})
                hdr = dicominfo(filelist{id});
                img = dicomread(hdr);
                
                
            elseif isNIFTI(filelist{id})
                nii = load_untouch_nii(filelist{id});
                img = nii.img;
                
            else
                errormsg = 'Unknown file type - DYNAMIC';
                return;
            end
            
            if ~isequal(size(img), size(TUMOR))
                errormsg = 'Unknown file type - DYNAMIC';
                return;
            end
            
            if i == 1
                
                DYNAMIC = img;
            else
                DYNAMIC(:,:,end+1:end+size(img,3)) = img;
            end
        end
    else
        for i = 1:size(LUT,1)
            id = LUT(i,1);
            if isDICOM(filelist{id})
                hdr = dicominfo(filelist{id});
                img = dicomread(hdr);
                
                
            elseif isNIFTI(filelist{id})
                nii = load_untouch_nii(filelist{id});
                img = nii.img;
                
            else
                errormsg = 'Unknown file type - DYNAMIC';
                return;
            end
            
            if ~isequal(size(img), size(TUMOR))
                errormsg = 'Unknown file type - DYNAMIC';
                return;
            end
            
            if i == 1
                
                DYNAMIC = img;
            else
                DYNAMIC(:,:,end+1:end+size(img,3)) = img;
            end
        end
    end
elseif filevolume == 3  
    % 2D slices
    for i = 1:size(LUT,1)
        for j = 1:size(LUT,2)
            id = LUT(i,j);
            if isDICOM(filelist{id})
                hdr = dicominfo(filelist{id});
                img = dicomread(hdr);
     
            elseif isNIFTI(filelist{id})
                nii = load_untouch_nii(filelist{id});
                img = nii.img;
                
            else
                errormsg = 'Unknown file type - DYNAMIC';
                return;
            end
            
            if ~isequal(size(img), size(TUMOR))
                errormsg = 'Unknown file type - DYNAMIC';
                return;
            end
            
            if i == 1 && j == 1
                
                DYNAMIC = img;
            else
                DYNAMIC(:,:,end+1:end+size(img,3)) = img;
            end
        end
    end
 
end
%% Resort DYNAMIC if fileorder is xytz

if ~strcmp(fileorder, 'xyzt')
    % reorder needed
    
    zslices = size(TUMOR,3);
    tslices = size(DYNAMIC,3)/zslices;
    
    for i = 1:tslices
        vect = [1:tslices:size(DYNAMIC,3)]+(i-1);
        
        curDYNAMIC = DYNAMIC(:,:,vect);
        
        newDYNAMIC(:,:,(zslices*(i-1)+[1:zslices]))=curDYNAMIC;
    end
end
        




















