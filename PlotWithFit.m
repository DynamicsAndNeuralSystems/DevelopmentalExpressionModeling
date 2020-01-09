function [xBinCenters,c10,cFree] = PlotWithFit(xData,yData,numBins,includeScatter,sizeIndicator)

f = figure('color','w');
hold('on')

%-------------------------------------------------------------------------------
% BINNING:
[xBinCenters,xThresholds,yMeans,yMedians] = makeQuantiles(xData,yData,numBins+1);
scatter(xBinCenters,yMeans,80*sizeIndicator,sizeIndicator,'filled')
plot(xBinCenters,yMeans,':k')
colormap([flipud(BF_getcmap('reds',9,0)); 0,0,0])

%-------------------------------------------------------------------------------
% Plot the data:
%-------------------------------------------------------------------------------
if includeScatter
    h_scatter = plot(xData,yData,'.k');
end
xlabel('Distance (m)')
ylabel('CGE')
f.Position = [1000        1012        1000         326];

%-------------------------------------------------------------------------------
% Exponential Fits:
%-------------------------------------------------------------------------------
xDataFitFull = xBinCenters'; yDataFitFull = yMeans';
% subset
includeMe = (sizeIndicator==1);
xDataFit = xDataFitFull(includeMe);
yDataFit = yDataFitFull(includeMe);

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

h_fit1 = plot(xDataFitFull,f_handle_exp10(xDataFitFull),'-r','LineWidth',3);
h_fit2 = plot(xDataFitFull,f_handle_expFree(xDataFitFull),'-g','LineWidth',3);

end
