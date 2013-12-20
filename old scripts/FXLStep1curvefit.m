%% Wrapper function for FXLStep1AIF.,

    function x = FXLStep1curvefit(xdata, i)

% Grab the relevant voxel i data

    Ct = xdata{1}.Ct;
    Ct = Ct(:,i); 
    xdata{1}.Ct = Ct;
