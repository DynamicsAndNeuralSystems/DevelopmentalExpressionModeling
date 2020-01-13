function [xBinCenters,xThresholds,yMeans,yMedians] = BF_PlotQuantiles(xData,yData,numThresholds,alsoScatter,makeNewFigure)
% Plots x-y scatter, but with mean of y plotted in quantiles of x
% Ben Fulcher
%-------------------------------------------------------------------------------

if nargin < 3 || isempty(numThresholds)
    numThresholds = 10;
end
if nargin < 4
    alsoScatter = false;
end
if nargin < 5
    makeNewFigure = false;
end

%-------------------------------------------------------------------------------
% Filter out NaNs:
goodBoth = (~isnan(xData) & ~isnan(yData));
if ~any(goodBoth)
    error('No good data');
elseif any(~goodBoth)
    xData = xData(goodBoth);
    yData = yData(goodBoth);
    fprintf(1,'Removed %u bad samples from x/y data\n',sum(~goodBoth));
end

xThresholds = arrayfun(@(x)quantile(xData,x),linspace(0,1,numThresholds));
xThresholds(end) = xThresholds(end) + eps; % make sure all data included in final bin
yMeans = arrayfun(@(x)mean(yData(xData>=xThresholds(x) & xData < xThresholds(x+1))),1:numThresholds-1);
yMedians = arrayfun(@(x)median(yData(xData>=xThresholds(x) & xData < xThresholds(x+1))),1:numThresholds-1);
yStds = arrayfun(@(x)std(yData(xData>=xThresholds(x) & xData < xThresholds(x+1))),1:numThresholds-1);

% ------------------------------------------------------------------------------
% Plot:
if makeNewFigure
    f = figure('color','w'); box('off');
end
hold('on');
theColor = [.59 .87 .82];
theStyle = '-';
theLineWidth = 1;

if alsoScatter
    plot(xData,yData,'o', 'MarkerSize', 10 ,'MarkerFaceColor', [.7 .7 .7], 'MarkerEdgeColor', [.6 .6 .6]);
end

xBinCenters = zeros(numThresholds-1,1);
for k = 1:numThresholds-1
    plot(xThresholds(k:k+1),ones(2,1)*yMeans(k),'LineStyle',theStyle,'LineWidth',theLineWidth,'Color', [0 0 0]);
    xBinCenters(k) = mean(xThresholds(k:k+1));
    plot(xBinCenters(k),yMeans(k),'o','MarkerSize',7,'LineStyle',theStyle,'LineWidth',theLineWidth,'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', [0 0 0]);
end

end
