%Approximate rule for integration. We are only given a certain number of
%sample points, so we need to determine the correct number of samples to
%fit in.

function S = sampleintegration(X,Y)

%X is the abscissae
%Y is the ordinate vector
S = 0;
if(numel(X) ~= numel(Y))
    error('X and Y values are not equal!');
elseif(numel(X) == 1)
    %1 value, integral is 0.
    S = 0;
elseif(numel(X) < 3)
    %2 values, use trapezoidal rule
    S = trapz(X, Y);
elseif(numel(X) < 4)
    %3 values, use Simpson's rule
    b = X(end);
    a = X(1);
    S = ((b-a)/6)*(Y(1)+4*Y(2)+Y(end));
elseif(numel(X) < 5)
    %4 values, use 3/8 Simpson's rule
    b = X(end);
    a = X(1);
    S = ((b-a)/8)*(Y(1) + 3*Y(2) + 3*Y(3) + Y(4));
elseif(numel(X) < 6)
    %5 values, use Boole's rule
    b = X(end);
    a = X(1);
    S = ((b-a)/90)*(7*Y(1) + 32*Y(2) + 12*Y(3) + 32*Y(4) + 7*Y(5));
    
else
    %Split the vector function into chunks for composite addition
    
    testX = X;
    testY = Y;
    
    total = numel(testX);
    fives = floor(total/5);
    remant = rem(total, 5);
 
    for k = 1:fives
        S = S+sampleintegration(testX(1:5), testY(1:5));
        testX(1:4) = [];
        testY(1:4) = [];
    end
   
    S = S + sampleintegration(testX, testY);
end
end
        
        
        
        
    
    
% else
%     %Split the vector function into chunks for composite simpson's and
%     %other stuff
%     
%     spareintervals = rem(numel(X) - 1, 3);
%     
%     if(spareintervals == 1)
%         %Find the smallest of front and back
%         
%         small = min(Y(1), Y(end));
%         
%         ind   = find(small == Y);
%         
%         if(ind == 1)
%         
%         S = sampleintegration(X(1:2), Y(1:2)) + sampleintegration(X(2:end), Y(2:end));
%         
%         else
%             
%               S = sampleintegration(X(end-1:end), Y(end-1:end)) + sampleintegration(X(1:end-1), Y(1:end-1));
%         end
%         
%     else
%         
%          S = sampleintegration(X(end-1:end), Y(end-1:end)) + sampleintegration(X(1:2), Y(1:2)) + sampleintegration(X(2:end-1), Y(2:end-1));
%     end
% end
% end
    