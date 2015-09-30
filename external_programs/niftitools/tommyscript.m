%Tommy script

%[FLASH] = batchextract('E:\data\exvivo++\', 'E:\data\Properexvivo++2');
%[RARE] = batchextract('E:\data\exvivo++\RARE\', 'E:\data\Properexvivo++2');
%[FLASHnot] = batchextract('E:\data\exvivonot++\', 'E:\data\Properexvivo++2');
%[RAREnot] = batchextract('E:\data\exvivonot++\RARE\', 'E:\data\Properexvivo++2');
%[returnRARE] = intensityscale('E:\temp\', 'E:\data\Properexvivo++2')

%a = mean(FLASH.means);
%b = mean(RARE.means);
%c = mean(FLASHnot.means);
%d = mean(RAREnot.means);



%d-b;
%c-a;
%b - a

%d - c
%[tester] = batchextract('E:\data\exvivo++\test\', 'C:\Documents and Settings\Thomas Ng\Desktop\lexvivo++2');
%[listing] = squishimage('E:\data\exvivo++\test\mwrxfl_xxxxstrip_Azmice2485_D71_04_hist')

%e = std(FLASH.means)
%f = std(RARE.means)
%g = std(FLASHnot.means)
%h = std(RAREnot.means)

%[numplaques] = blackdot('E:\data\exvivo++\mwrxfl_xxxxstrip_Azmice1570_Dc1_04roiextracted', 'E:\data\exvivo++\Properexvivo++2');

%Perform unsharp filtering
 
%realpath1 = strcat('E:\temp\', '*scaled.hdr')
% listoffiles1 = ls(realpath1)
% sizer1 = size(listoffiles1)
% length2 = sizer1(1)
 
 
 
 % for i=1:length2
  %   currentpath = strcat('E:\temp\', listoffiles1(i,:))
   %  a = load_nii(currentpath)     
  
    % sharpFilter = fspecial('unsharp', 0.2);
     
    % sharp = imfilter(a.img, sharpFilter, 'conv');
      
    % a.img = sharp;
     
    %filterlist = returnRARE.filtername;
     
     %save_nii(a, filterlist(2*i,:));
 %end

cd /home/tommy/projects/AZ/matlab/

%pathscript

load('exvivoplaque')
template = load_nii('/media/windows2/data/Properexvivo++2');
templatesize = find(template.img > 0);
templatesize = size(templatesize);
templatevoxels = templatesize(1);

%in mm
voxdim = [0.05 0.05 0.05];

done = extracttables(exvivoplaque, 3, ['/home/tommy/projects/AZ/data/' ...
                    'histograms/'], templatevoxels, voxdim);

%Draw Histogram

patientno = size(done);
patientno = patientno(2);

for i = 1:patientno
  voxsize = done(i).plaquemmvolume;
  histogram(i,:) = voxsize.*done(i).histogram;
  percent(i) = done(i).percentplaque;
  slabs186(i,:) = voxsize.*done(i).slabhist186;
  slabs228(i,:) = voxsize.*done(i).slabhist228;
  slabs450(i,:) = voxsize.*done(i).slabhist450;
  
end
subplot(1,2,1);
h = bar(histogram, 'stack');

display = '';
for i = 1:patientno
  set(get(h(i), 'DisplayName'), 'DisplayName', done(i).name);
  
  
display = strvcat(display,strcat(num2str(i) , '--', done(i).name));  
  
  
end

slabdisplay = strvcat('1 -- 152 to 186', '2 -- 187 to 228', ['3 -- 228 to ' ...
                    'end']);

title('Plaque Volume of whole ROI and proportion of intensity')
ylabel('Volume of Plaque (mm^3)')
legend(' p = [0.01, 0.005) - soft', 'p = [0.005, 0.001) - medium ', ['c = [0.001,' ...
                    'infinity) - large'])
xlabel(display)

subplot(1,2,2);
f = bar(percent);
title('Plaque Volume of whole ROI as Percentage')
ylabel('Percentage of total ROI volume')
xlabel(display)

figure;

subplot(2,3,1);
h = bar(slabs186, 'stack');
title('Plaque Volume Slices 152 to 186')
ylabel('Volume of Plaque (mm^3)')
axis([0,4,0,(6.5)*10^4]);
xlabel(display)

subplot(2,3,2);
i = bar(slabs228, 'stack');
title('Plaque Volume Slices 187 to 228')
ylabel('Volume of Plaque (mm^3)')
axis([0,4,0,(6.5)*10^4]);
xlabel(display)
legend(' p = [0.01, 0.005) - soft', 'p = [0.005, 0.001) - medium ', ['c = [0.001,' ...
                    'infinity) - large'])

subplot(2,3,3);
j = bar(slabs450, 'stack');

title('Plaque Volume Slices 228 to end')
ylabel('Volume of Plaque (mm^3)')
xlabel(display)
axis([0,4,0,(6.5)*10^4]);


subplot(2,3,4);
k = bar([slabs186(1,:);slabs228(1,:);slabs450(1,:)], 'stack');
g = done(1).name;
title(strcat('Plaque Volume for Mouse', g));

ylabel('Volume of Plaque (mm^3)')
axis([0,4,0,(6.5)*10^4]);
xlabel(slabdisplay)


subplot(2,3,5);
k = bar([slabs186(2,:);slabs228(2,:);slabs450(2,:)], 'stack');
g = done(2).name;
title(strcat('Plaque Volume for Mouse', g));

ylabel('Volume of Plaque (mm^3)')
axis([0,4,0,(6.5)*10^4]);
xlabel(slabdisplay)

subplot(2,3,6);
k = bar([slabs186(3,:);slabs228(3,:);slabs450(3,:)], 'stack');
g = done(3).name;
title(strcat('Plaque Volume for Mouse', g));

ylabel('Volume of Plaque (mm^3)')
axis([0,4,0,(6.5)*10^4]);
xlabel(slabdisplay)

for i = 1:patientno
  
  done(i).percent186 = done(i).slabhist186 / sum(done(i).slabhist186)
  done(i).percent228 = done(i).slabhist228 / sum(done(i).slabhist228);
  done(i).percent450 = done(i).slabhist450 / sum(done(i).slabhist450);
  
  done(i).percenttotal186 = done(i).slabhist186 / done(i).plaquevoxelvolume;
  done(i).percenttotal228 = done(i).slabhist228 / done(i).plaquevoxelvolume;
  done(i).percenttotal450 = done(i).slabhist450 / ...
      done(i).plaquevoxelvolume;
  
  
  
end