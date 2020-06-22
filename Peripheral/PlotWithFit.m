function [xBinCenters,c10,cFree] = PlotWithFit(xData,yData,numBins,includeScatter,sizeIndicator,newFig)
% Plot x-y data across a set of equiprobable bins (and fitted exponential)
%-------------------------------------------------------------------------------
if nargin < 6
    newFig = false;
end

if newFig
    f = figure('color','w');
    hold('on')
end

%-------------------------------------------------------------------------------
% BINNING:
[xBinCenters,xThresholds,yMeans,yMedians] = makeQuantiles(xData,yData,numBins+1);
scatter(xBinCenters,yMeans,40*sizeIndicator,sizeIndicator)
plot(xBinCenters,yMeans,':k')
colormap([flipud(BF_getcmap('reds',9,0)); 0,0,0])

% Bin extents:
hold('on')
for k = 1:numBins
    plot(xThresholds(k:k+1),ones(2,1)*yMeans(k),'LineStyle','-',...
        'LineWidth',1,'Color','k')
    % plot(mean(xThresholds(k:k+1)),yBinMeans(k),'o','MarkerSize',5,'LineStyle',theStyle,...
    %     'LineWidth',theLineWidth,'Color',theColor)
end


%-------------------------------------------------------------------------------
% Plot the data:
%-------------------------------------------------------------------------------
if includeScatter
    h_scatter = plot(xData,yData,'.k');
end
xlabel('Distance (mm)')
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

doExp10 = false;

% A--exp10
if doExp10
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint',1/5);
    f = fittype('exp(-x*n)','options',s);
    c10 = fit(xDataFit,yDataFit,f);
    f_handle_exp10 = @(x) exp(-x*c10.n);
end

% B--expFree
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1,1/5,0]);
f = fittype('A*exp(-x*n) + B','options',s);
cFree = fit(xDataFit,yDataFit,f);
f_handle_expFree = @(x) cFree.A*exp(-x*cFree.n) + cFree.B;

if doExp10
    h_fit1 = plot(xDataFitFull,f_handle_exp10(xDataFitFull),'-r','LineWidth',3);
end
h_fit2 = plot(xDataFitFull,f_handle_expFree(xDataFitFull),'-g','LineWidth',3);


end
