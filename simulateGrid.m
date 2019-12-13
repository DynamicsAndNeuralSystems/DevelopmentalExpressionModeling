% Make a grid in 2d space
whatGradients = 'ExpDecaySingle'; % 'linear' 'poly', 'Gaussian', 'GaussianFixedScale'
numGradients = 200;

% Make 2-d spatial grid in [X,Y]:
extentSim = 30;
extentZoom = 10;
resolution = 20;
[coOrds,X,Y] = MakeGrid(extentSim,resolution);
xRange = [min(X(:)),max(X(:))];
yRange = [min(Y(:)),max(Y(:))];
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
        X0 = rand*(xRange(2)-xRange(1))+xRange(1);
        Y0 = rand*(yRange(2)-yRange(1))+yRange(1);
        rCoeff = rand(2,1)-0.5;
        g = rCoeff(1)*(X-X0) + rCoeff(2)*(Y-Y0);
    case 'poly'
        pOrderX = randi(1);
        pOrderY = randi(2);
        X0 = rand*(xRange(2)-xRange(1))+xRange(1);
        Y0 = rand*(yRange(2)-yRange(1))+yRange(1);
        rCoeff = rand(2,1)-0.5;
        g = rCoeff(1)*(X-X0).^pOrderX + rCoeff(2)*(Y-Y0).^pOrderY;
    case 'Gaussian'
        A = rand();
        X0 = rand*(xRange(2)-xRange(1))+xRange(1);
        Y0 = rand*(yRange(2)-yRange(1))+yRange(1);
        sigma2X = rand*(xRange(2)-xRange(1))/2;
        sigma2Y = rand*(yRange(2)-yRange(1))/2;
        g = A*exp(-((X-X0).^2./sigma2X + (Y-Y0).^2./sigma2Y));
    case 'GaussianFixedScale'
        A = 1;
        scaleProp = 1/4;
        X0 = rand*(xRange(2)-xRange(1))+xRange(1);
        Y0 = rand*(yRange(2)-yRange(1))+yRange(1);
        sigma2X = (xRange(2)-xRange(1))*scaleProp;
        sigma2Y = (yRange(2)-yRange(1))*scaleProp;
        g = A*exp(-(((X-X0)/sigma2X).^2 + ((Y-Y0)/sigma2Y).^2));
    case 'ExpDecaySingle'
        d0 = rand*20;
        X0 = rand*(xRange(2)-xRange(1))+xRange(1);
        Y0 = rand*(yRange(2)-yRange(1))+yRange(1);
        D = (X-X0).^2 + (Y-Y0).^2;
        g = exp(-D/d0) + 0.2*randn(size(D));
    end
    gStretch = g(:);
    expData(:,i) = gStretch;
end

%-------------------------------------------------------------------------------
% Now we zoom in
doZoom = false;
if doZoom
    extentRange = [extentSim/2-extentZoom/2,extentSim/2+extentZoom/2];
    isInRange = @(x) x>=extentRange(1) & x<=extentRange(2);
    keepMe = isInRange(coOrds(:,1)) & isInRange(coOrds(:,2));
    expDataZoom = expData(keepMe,:);
    coOrdsZoom = coOrds(keepMe,:);
    XKeep = isInRange(X(1,:));
    YKeep = isInRange(Y(:,1));
    XZoomed = X(YKeep,XKeep);
    YZoomed = Y(YKeep,XKeep);
else
    expDataZoom = expData;
    coOrdsZoom = coOrds;
    XZoomed = X;
    YZoomed = Y;
end

% Compute pairwise similarity of gradients:
expDataZ = zscore(expDataZoom);
cgeVect = 1 - pdist(expData,'corr');

% Compute pairwise distance matrix:
dVectZoom = pdist(dScale*coOrdsZoom,'Euclidean');
dMatZoom = squareform(dVectZoom);

%-------------------------------------------------------------------------------
% Plot some individual gradients:
f = figure('color','w');
for i = 1:12
    subplot(3,4,i)
    f.Position = [1000        1012        1000         326];
    % Plot as a grid:
    imagesc(reshape(expData(:,i),size(XZoomed)))
    axis('square')
end
colormap(gray)

%-------------------------------------------------------------------------------
% Plot against each other:
f = figure('color','w');
f.Position = [1000        1012        1000         326];
numBins = 50;
plot(dVectZoom,cgeVect,'.k')
hold('on')
[xBinCenters,xThresholds,yMeans,yMedians] = BF_PlotQuantiles(dVectZoom(:),cgeVect(:),numBins);
xlabel('Distance (m)')
ylabel('CGE')
title(sprintf('%u superimposed %s gradients',numGradients,whatGradients))

% FIT EXPONENTIAL:
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[0.3,1/5,0]);
f = fittype('A*exp(-x*n) + B','options',s);
xData = xBinCenters; yData = yMeans';
% xData = dVectZoom; yData = cgeVect;
[c, Stats] = fit(xData,yData,f);
f_handle = @(x) c.A.*exp(-x*c.n) + c.B;

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
