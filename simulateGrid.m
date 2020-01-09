% Make a grid in 2d space
whatGradients = 'spatialLag'; % 'linear' 'poly', 'Gaussian', 'GaussianFixedScale', 'ExpDecaySingle'
numGradients = 1000;

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
        d0 = 5; %(max(coOrds(:,1))-min(coOrds(:,1)));
        rho = 0.25;
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
% A--exp10
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',1/5);
f = fittype('exp(-x*n)','options',s);
xData = xBinCenters; yData = yMeans';
% xData = dVectZoom; yData = cgeVect;
c10 = fit(xData,yData,f);
f_handle_exp10 = @(x) exp(-x*c10.n);

% B--expFree
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1,1/5,0]);
f = fittype('A*exp(-x*n) + B','options',s);
xData = xBinCenters; yData = yMeans';
% xData = dVectZoom; yData = cgeVect;
cFree = fit(xData,yData,f);
f_handle_expFree = @(x) cFree.A*exp(-x*cFree.n) + cFree.B;

hold('on')
plot(xData,f_handle_exp10(xData),'-r','LineWidth',3)
plot(xData,f_handle_expFree(xData),'-g','LineWidth',3)
title(sprintf('%u superimposed gradients: n = %g',numGradients,1/c.n))

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
scalingFactor = [1,1.5,2,2.5,3];
numThresholds = 20; % for histogram/fitting

% Fitting:
paramStruct = cell(length(scalingFactor),1);
for i = 1:length(scalingFactor)
    [xBinCenters,xThresholds,yMeans,yMedians] = makeQuantiles(scalingFactor(i)*dVectZoom(:),cgeVect(:),numThresholds);
    xData = xBinCenters'; yData = yMeans';
    [Exp_3free_fun,Stats,paramStruct{i}] = GiveMeFit(xData,yData,'exp',true);
end

% Plotting:
f = figure('color','w');
colors = BF_getcmap('spectral',5,1);
hold('on');
for i = 1:length(scalingFactor)
    [xBinCenters,xThresholds,yMeans,yMedians] = makeQuantiles(scalingFactor(i)*dVect,cgeVect,numThresholds);
    plot(xBinCenters,yMeans,'o','color',colors{i})
    plot(xBinCenters,Exp_3free_fun(xBinCenters,paramStruct{i}),'-','color',colors{i},'LineWidth',2)
end

% Plot parametric variation:
f = figure('color','w');
subplot(3,1,1);
lambdaEst = arrayfun(@(x)1/paramStruct{x}.n,1:length(paramStruct))
plot(scalingFactor,lambdaEst,'o-k')
title('Spatial scale')
subplot(3,1,2);
strengthEst = arrayfun(@(x)paramStruct{x}.A,1:length(paramStruct))
plot(scalingFactor,strengthEst,'o-k')
title('Strength')
subplot(3,1,3);
offSetEst = arrayfun(@(x)paramStruct{x}.B,1:length(paramStruct))
plot(scalingFactor,offSetEst,'o-k')
title('Offset')
