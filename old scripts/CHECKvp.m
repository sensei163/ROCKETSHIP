%% Check the fit

function check = CHECKvp(x, xdata)

r = round(1 + (size(x,1)-1).*rand(min(9,size(x,1)),1));
a = figure;
%r = [688 1557 1664 1999];
for j = 1:numel(r)
    i = r(j);
    timer = xdata{1}.timer;
    Ct = xdata{1}.Ct;
    Ct = Ct(:,i);
    
    %%FIX below
    y(j,:) = FXLStep1AIFcfit(x(i,1), x(i,2), xdata{1}.Cp, xdata{1}.timer);
    
    
    figure(a)
    subplot(3,3,j), plot(timer, Ct, 'x'), hold on, plot(timer, y(j,:), 'r'), title(['Vo:' num2str(i), 'Ktve:' num2str(x(i,1:2)) 'R2:' num2str(x(i,4))])
    hold off
end

check = 1;