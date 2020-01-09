function [c10,cFree] = PlotWithFit(xData,yData,numBins,includeScatter)

f = figure('color','w');
hold('on')

%-------------------------------------------------------------------------------
% BINNING:
[xBinCenters,xThresholds,yMeans,yMedians] = BF_PlotQuantiles(xData(:),yData(:),numBins);

%-------------------------------------------------------------------------------
% Plot the data:
%-------------------------------------------------------------------------------
h_scatter = plot(xData,yData,'.k');
xlabel('Distance (m)')
ylabel('CGE')
f.Position = [1000        1012        1000         326];

%-------------------------------------------------------------------------------
% Exponential Fits:
%-------------------------------------------------------------------------------
xDataFit = xBinCenters; yDataFit = yMeans';
% xData = dVectZoom; yData = cgeVect;

% A--exp10
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',1/5);
f = fittype('exp(-x*n)','options',s);
c10 = fit(xDataFit,yDataFit,f);
f_handle_exp10 = @(x) exp(-x*c10.n);

% B--expFree
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1,1/5,0]);
f = fittype('A*exp(-x*n) + B','options',s);
cFree = fit(xDataFit,yDataFit,f);
f_handle_expFree = @(x) cFree.A*exp(-x*cFree.n) + cFree.B;

h_fit1 = plot(xDataFit,f_handle_exp10(xDataFit),'-r','LineWidth',3);
h_fit2 = plot(xDataFit,f_handle_expFree(xDataFit),'-g','LineWidth',3);

end
