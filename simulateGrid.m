% Make a grid in 2d space
whatGradients = 'spatialLag'; % 'linear' 'poly', 'Gaussian', 'GaussianFixedScale'
numGradients = 100;

% Make 2-d spatial grid in [X,Y]:
extentSim = 50;
extentZoom = 20;
resolution = 25;
[coOrds,X,Y] = MakeGrid(extentSim,resolution);
numAreas = length(coOrds);

% Distance (m) corresponding to a unit change in each axis
dScale = 1;

% Compute pairwise distance matrix:
dVect = pdist(dScale*coOrds,'Euclidean');
dMat = squareform(dVect);
%-------------------------------------------------------------------------------

% Compute gradients through the space:
expData = zeros(numAreas,numGradients);
for i = 1:numGradients
    switch whatGradients
    case 'spatialLag'
        d0 = 2; %(max(coOrds(:,1))-min(coOrds(:,1)));
        rho = 0.5;
        g = GenerateSpatialLagMap(dMat,d0,rho);
    case 'linear'
        X0 = rand*(xRange(end)-xRange(1))+xRange(1);
        Y0 = rand*(yRange(end)-yRange(1))+yRange(1);
        rCoeff = rand(2,1)-0.5;
        g = rCoeff(1)*(X-X0) + rCoeff(2)*(Y-Y0);
    case 'poly'
        pOrderX = randi(1);
        pOrderY = randi(2);
        X0 = rand*(xRange(end)-xRange(1))+xRange(1);
        Y0 = rand*(yRange(end)-yRange(1))+yRange(1);
        rCoeff = rand(2,1)-0.5;
        g = rCoeff(1)*(X-X0).^pOrderX + rCoeff(2)*(Y-Y0).^pOrderY;
    case 'Gaussian'
        A = rand();
        X0 = rand*(xRange(end)-xRange(1))+xRange(1);
        Y0 = rand*(yRange(end)-yRange(1))+yRange(1);
        sigma2X = rand*(xRange(end)-xRange(1))/2;
        sigma2Y = rand*(yRange(end)-yRange(1))/2;
        g = A*exp(-((X-X0).^2./sigma2X + (Y-Y0).^2./sigma2Y));
    case 'GaussianFixedScale'
        A = 1;
        X0 = rand*(xRange(end)-xRange(1))+xRange(1);
        Y0 = rand*(yRange(end)-yRange(1))+yRange(1);
        sigma2X = (xRange(end)-xRange(1))/5;
        sigma2Y = (yRange(end)-yRange(1))/5;
        g = A*exp(-((X-X0).^2./sigma2X + (Y-Y0).^2./sigma2Y));
    end
    gStretch = g(:);
    expData(:,i) = gStretch;
end
expDataZ = zscore(expData);

% Compute pairwise similarity of gradients:
cgeVect = 1 - pdist(expDataZ,'corr');

%-------------------------------------------------------------------------------
% Plot some individual gradients:
f = figure('color','w');
for i = 1:12
    subplot(3,4,i)
    f.Position = [1000        1012        1000         326];
    % Plot as a grid:
    imagesc(reshape(expDataZ(:,i),size(X)))
    axis('square')
end


%-------------------------------------------------------------------------------
% Plot against each other:
f = figure('color','w');
f.Position = [1000        1012        1000         326];
numBins = 25;
plot(dVect,cgeVect,'.k')
hold('on')
[xBinCenters,xThresholds,yMeans,yMedians] = BF_PlotQuantiles(dVect(:),cgeVect(:),numBins);
xlabel('Distance (m)')
ylabel('CGE')
title(sprintf('%u superimposed %s gradients',numGradients,whatGradients))

% FIT EXPONENTIAL:
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[rho,d0,0]);
f = fittype('A*exp(-x/n) + B','options',s);
% xData = xBinCenters; yData = yMeans';
xData = dVect(:); yData = cgeVect(:);
[c, Stats] = fit(xData,yData,f);
f_handle = @(x) c.A.*exp(-x/c.n) + c.B;

% [f_handle,Stats,c] = GiveMeFit(xBinCenters,yMeans,'exp',false);

% Plot as a grid:
% subplot(1,2,1)

% [xBinCenters,xThresholds,yMeans,yMedians] = BF_PlotQuantiles(xData,yData,numThresholds,alsoScatter,makeNewFigure)
% surf(X,Y,reshape(sum(expDataZ,2),size(X)))
% xlabel('x')
% ylabel('y')
% zlabel('z')
%-------------------------------------------------------------------------------


%-------------------------------------------------------------------------------
% Now stretch it out uniformly:
numThresholds = 10;
f = figure('color','w');
hold('on');
scalingFactor = [1,1.2,1.5,2,3];
colors = BF_getcmap('spectral',5,1);
for i = 1:length(scalingFactor)
    [xBinCenters,xThresholds,yMeans,yMedians] = makeQuantiles(scalingFactor(i)*dVect,cgeVect,numThresholds);
    plot(xBinCenters,yMeans,'o-','color',colors{i})
end
