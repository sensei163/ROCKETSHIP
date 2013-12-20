%% Helper Function for nonlinear curvefit to FXLAIF model, no vp
function GG = FXRAIFhelper_A(xdata, voxel)
warning off
    
R1t = xdata{1}.R1t;
i  = voxel;
R1t = R1t(:,i);

Cp = xdata{1}.Cp;

R1o = xdata{1}.R10;
R1i = xdata{1}.R1i;

fw  = xdata{1}.fw;
r1  = xdata{1}.r1;

% Use Curvefitting tool box

options = fitoptions('Method', 'NonlinearLeastSquares', 'Algorithm', 'Levenberg-Marquardt', 'MaxIter', 150, 'MaxFunEvals', 150, 'TolFun', 0.001, 'TolX', 0.001, 'Display', 'off', 'Lower',[0 0 0], 'Upper', [5 1 5], 'StartPoint', [0.5 0.5 0.5]);


ft = fittype('FXRAIFtest(Ktrans, ve, tau1, Cp, T1)', 'independent', {'T1', 'Cp'},'coefficients',{'Ktrans', 've', 'tau1'});



  timer = xdata{1}.timer';


[f gof] = fit([timer, xdata{1}.Cp],R1t,ft, options)

toc
%{
%Calculate the R2 fit
GG(1) = f.Ktrans;
GG(2) = f.ve;
GG(3) = f.tau1;
GG(end+1) = gof.rsquare;
GG(end+1) = gof.adjrsquare;
%}
GG = 1;

GG = [f.Ktrans f.ve f.tau1 gof.rsquare gof.adjrsquare];

