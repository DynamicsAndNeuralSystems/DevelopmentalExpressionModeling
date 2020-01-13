%-------------------------------------------------------------------------------
% Parameters:
%-------------------------------------------------------------------------------
% Parameters of the 2D spatial coordinate grid:
extentSim = 50; % physical extent of the
resolution = 40; % number of points to split each dimension into
dScale = 1; % distance (m) corresponding to a unit change in each axis

% Details of the expression maps to generate (what parameters of what model):
whatGradients = 'spatialLag'; % 'linear' 'poly', 'Gaussian', 'GaussianFixedScale', 'ExpDecaySingle'
numGradients = 100;

% Spatial binning for plotting and analysis:
numBins = 25;

% Zoom in for analysis:
doZoom = false;
extentZoom = 10;

%-------------------------------------------------------------------------------
% Generate the 2-d spatial grid in [X,Y] (on which to work):
[coOrds,X,Y] = MakeGrid(extentSim,resolution);
numAreas = length(coOrds);

% Compute the pairwise distance matrix in the coordinate space:
dVect = pdist(dScale*coOrds,'Euclidean');
dMat = squareform(dVect);

%-------------------------------------------------------------------------------
% Understand the spatial sampling of points across distance bins:
propRegionsRepresented = AnalyseSpatialSampling(dVect,numAreas,numBins);

%-------------------------------------------------------------------------------
% Generate a bunch of spatial gene-expression gradients:
expData = GenerateEnsemble(whatGradients,numAreas,numGradients,dMat);

%-------------------------------------------------------------------------------
% Zoom in on a subset (e.g., for finite-size effects):
if doZoom
    [expDataZoom,coOrdsZoom,XZoomed,YZoomed] = ZoomIn(extentSim,extentZoom,coOrds,expData,X,Y);
else
    expDataZoom = expData;
    coOrdsZoom = coOrds;
    XZoomed = X;
    YZoomed = Y;
end

% Compute the pairwise distance matrix
% (in matrix and vector form):
dVectZoom = pdist(dScale*coOrdsZoom,'Euclidean');
dMatZoom = squareform(dVectZoom);

%-------------------------------------------------------------------------------
% Compute pairwise similarity of gradients as CGE:
% Normalize each gradient:
normalizeHow = 'zscore'; % 'zscore','mixedSigmoid','subtractMean'
expDataNorm = BF_NormalizeMatrix(expDataZoom,normalizeHow);
% expDataNorm = expDataZoom;
% for i = 1:size(expDataZoom,2)
%     expDataNorm(:,i) = expDataNorm(:,i)+randn(1);
% end

cgeVectNorm = 1 - pdist(expDataNorm,'corr');
cgeVect = 1 - pdist(expDataZoom,'corr');

%-------------------------------------------------------------------------------
% Plot some individual gradients as color in 2d:
PlotSomeIndividualGradients(expDataZoom,size(XZoomed));

%-------------------------------------------------------------------------------
% Plot:
includeScatter = false;

% Unnormalized data:
% [binCenters,c10,cFree] = PlotWithFit(dVectZoom,cgeVect,numBins,includeScatter,propRegionsRepresented)
% title(sprintf('%u superimposed %s gradients: n = %g',numGradients,whatGradients,1/cFree.n))

% Normalized data:
[binCenters,c10,cFree] = PlotWithFit(dVectZoom,cgeVectNorm,numBins,includeScatter,propRegionsRepresented)
title(sprintf('%u superimposed %s gradients: n = %g',numGradients,whatGradients,1/cFree.n))
% bar(binCenters,propRegionsRepresented)

% Together:
% f = figure('color','w'); hold('on')
% BF_PlotQuantiles(dVectZoom,cgeVect,numBins);
% BF_PlotQuantiles(dVectZoom,cgeVectNorm,numBins);
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
