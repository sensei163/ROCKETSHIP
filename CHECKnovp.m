%% Check the fit

function check = CHECKnovp(x, xdata)

r = round(1 + (size(x,1)-1).*rand(min(25,size(x,1)),1));
a = figure;

for j = 1:numel(r)
    i = r(j);
    timer = xdata{1}.timer;
    Ct = xdata{1}.Ct;
    Ct = Ct(:,i);
    
    y(j,:) = FXLStep1AIFcfit(x(i,1), x(i,2), xdata{1}.Cp, xdata{1}.timer);
    
    
    figure(a)
    subplot(5,5,j), plot(timer, Ct, 'x'), hold on, plot(timer, y(j,:), 'r'), title(['Voxel: ' num2str(i), 'Ktrans and ve: ' num2str(x(i,1:2)) 'R^2 fit: ' num2str(x(i,3))])
    hold off
end

check = 1;