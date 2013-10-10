%% Check the fit

function check = CHECKFXR(x, xdata)

r = round(1 + (size(x,1)-1).*rand(min(16,size(x,1)),1));

%r = [688 1557 1664 1999]
a = figure;

for j = 1:numel(r)
    i = r(j);
    timer = xdata{1}.timer;
    Rt1 = xdata{1}.R1tTOI;
    Rt1 = Rt1(:,i);
    
    y(j,:) =  FXRAIF(x(i,:), xdata); %FXLStep1AIFcfit(x(i,1), x(i,2), xdata{1}.Cp, xdata{1}.timer);
    
    
    figure(a)
    subplot(4,4,j), plot(timer, Rt1, 'x'), hold on, plot(timer, y(j,:), 'r'), title(['Voxel: ' num2str(i), 'K_vt: ' num2str(x(i,1:3)) 'R^2 fit: ' num2str(x(i,end))])
    hold off
end

check = 1;