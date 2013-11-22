load('C:\Users\sbarnes\Documents\data\6 DCE Stroke\sb01_06nov13.mH1\20131106_000_sb01_06nov13_dynamic_01_s_AIF_with_vpFIT_ROI.mataif_vp_FIT_voxels.mat')

voi = 99;
% voi = 1000;
% voi = 313;
% voi = 558;
curve = xdata{1,1}.Ct(:,voi);
moving_curve = smooth(curve,8,'moving');
lowess_curve = smooth(curve,0.01,'rlowess');

figure(1);
foo = FXLStep1AIFhelper(xdata,voi,0)
bar = FXLStep1AIF(foo(1:3),xdata);
plot(1:800,lowess_curve,1:800,bar(1,:));

figure(2);
foo_lowess = FXLStep1AIFhelper(xdata,voi,2)
bar_lowess = FXLStep1AIF(foo_lowess(1:3),xdata);
plot(1:800,lowess_curve,1:800,bar_lowess(1,:));

figure(3);
foo_moving = FXLStep1AIFhelper(xdata,voi,1)
bar_moving = FXLStep1AIF(foo_moving(1:3),xdata);
plot(1:800,moving_curve,1:800,bar_moving(1,:));

figure(4);
% plot(1:800,lowess_curve,1:800,bar(1,:),1:800,bar_lowess(1,:),1:800,bar_moving(1,:));
plot(1:800,bar(1,:),1:800,bar_lowess(1,:),1:800,bar_moving(1,:));
legend('normal','lowess','moving');



figure(5);
bar_mine = FXLStep1AIF(foo_mine(1:3),xdata);
plot(1:800,curve,1:800,bar_mine(1,:));



figure(3);
foo = NaN(size(currentimg));
% foo = zeros(size(currentimg));
foo(tumind) = xdata{1,1}.Ct(200,:);
imshow(foo.*20000)
h = fspecial('average', [2 2]);
bar = filter2(h, foo);
imshow(bar.*20000)