function x = cleanupr(x, rfit)

ind = find(x(:,end) < rfit);

x(ind,:) = [];
