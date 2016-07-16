nii = load_untouch_nii('Dynamic_t1w.nii');
img = nii.img;
img_org = nii.img;
figure, imshow3D(squeeze(img_org), [0 1200]);

fixed = img(:,:,:,1);
[optimizer,metric] = imregconfig('multimodal');
p = ProgressBar(size(img,4)-2);
for n=2:size(img,4)
    moving = img(:,:,:,n);
    movingRegistered = imregister(moving, fixed, 'affine', optimizer, metric);
    img(:,:,:,n) = movingRegistered;
    p.progress;
end
p.stop;



figure, imshow3D(squeeze(img), [0 1200]);
%nii.img = img;
%save_untouch_nii(nii,'registered.nii.gz')

% figure, imshowpair(moving, fixed, 'montage')
% title('Unregistered')
% 
% figure, imshowpair(moving, fixed)
% title('Unregistered')
% 
% [optimizer,metric] = imregconfig('multimodal');
% movingRegisteredDefault = imregister(moving, fixed, 'affine', optimizer, metric);
% 
% figure, imshowpair(movingRegisteredDefault, fixed)
% title('A: Default registration')