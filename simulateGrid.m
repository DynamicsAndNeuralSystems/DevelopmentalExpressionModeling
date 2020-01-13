%-------------------------------------------------------------------------------
% Parameters:
%-------------------------------------------------------------------------------
% Parameters of the 2D spatial coordinate grid:
extentSim = 50; % physical extent of the
resolution = 40; % number of points to split each dimension into
dScale = 1; % distance (m) corresponding to a unit change in each axis

% Details of the expression maps to generate (what parameters of what model):
whatGradients = 'spatialLag'; % 'linear' 'poly', 'Gaussian', 'GaussianFixedScale', 'ExpDecaySingle'
ensembleParams = struct();
ensembleParams.rho = 0.5; % relative strength of spatial autocorrelation (relative to noise)
ensembleParams.d0 = 1; % spatial scale
numGradients = 100;
% Normalization of each gradient:
normalizeHow = 'zscore'; % 'zscore','mixedSigmoid','subtractMean'

% Spatial binning for plotting and analysis:
numBins = 25;
includeScatter = false;

% Time Variation:
% --A combination of isotropic spatial scaling, and 'zooming out'--
scalingFactors = [1,1.5,2,2.5,3];
zoomFactors = [1,1,1,1,1];
assert(length(scalingFactors)==length(zoomFactors));
assert(all(zoomFactors >= 1));
numTimePoints = length(zoomFactors);
timeVector = 1:numTimePoints;

%-------------------------------------------------------------------------------
% Generate the 2-d spatial grid in [X,Y] (on which to work):
[coOrds,X,Y] = MakeGrid(extentSim,resolution);
numAreas = length(coOrds);

% Compute the pairwise distance matrix in the coordinate space:
dVect = pdist(dScale*coOrds,'Euclidean');
dMat = squareform(dVect);

%-------------------------------------------------------------------------------
% Understand the spatial sampling of points across distance bins:
propRegionsRepresented = AnalyseSpatialSampling(dMat,numAreas,numBins);

%-------------------------------------------------------------------------------
% Generate a bunch of spatial gene-expression gradients:
expData = GenerateEnsemble(whatGradients,numAreas,numGradients,dMat,ensembleParams);

%-------------------------------------------------------------------------------
% Plot some individual gradients as color in 2d:
PlotSomeIndividualGradients(expData,size(X));

%-------------------------------------------------------------------------------
% Now we go through 'time':
f = figure('color','w');
colors = BF_getcmap('spectral',numTimePoints,1);
hold('on')
paramStruct = cell(numTimePoints,1);
for t = timeVector
    % Scale the original map by this much (isotropic stretch):
    scalingFactor = scalingFactors(t);
    % Zoom out from the original map (with the same rules) by this much:
    zoomFactor = zoomFactors(t);

    %-------------------------------------------------------------------------------
    % Zoom in on a subset (e.g., for finite-size effects):
    if zoomFactor > 1
        [expDataT,coOrdsT,XT,YT] = ZoomIn(extentSim,extentZoom,coOrds,expData,X,Y);
    else
        % No zoom at this time point (just copy original variables)
        expDataT = expData;
        coOrdsT = coOrds;
        XT = X;
        YT = Y;
    end

    %-------------------------------------------------------------------------------
    % Rescale:
    XT = XT*scalingFactor;
    YT = YT*scalingFactor;
    coOrdsT = coOrdsT*scalingFactor;

    %-------------------------------------------------------------------------------
    % Recompute the pairwise distance matrix for this time point:
    % (in matrix and vector form):
    dVectT = pdist(scalingFactor*coOrdsT,'Euclidean');
    dMatT = squareform(dVectT);

    %-------------------------------------------------------------------------------
    % Normalize each gradient:
    expDataTNorm = BF_NormalizeMatrix(expDataT,normalizeHow);

    %-------------------------------------------------------------------------------
    % Compute pairwise similarity of gradients as CGE:
    cgeVectNorm = 1 - pdist(expDataTNorm,'corr');
    % (unnormalized:)
    cgeVect = 1 - pdist(expDataT,'corr');

    %-------------------------------------------------------------------------------
    % Fit:
    [xBinCenters,xThresholds,yMeans,yMedians] = makeQuantiles(dVectT,cgeVectNorm,numBins+1);
    [Exp_3free_fun,Stats,paramStruct{t}] = GiveMeFit(xBinCenters',yMeans','exp',true);

    %-------------------------------------------------------------------------------
    % Plot:
    % Unnormalized data:
    % [binCenters,c10,cFree] = PlotWithFit(dVectT,cgeVect,numBins,includeScatter,propRegionsRepresented)
    % title(sprintf('%u superimposed %s gradients: n = %g',numGradients,whatGradients,1/cFree.n))

    % Normalized data:
    subplot(numTimePoints,1,t);
    hold('on')
    [binCenters,c10,cFree] = PlotWithFit(dVectT,cgeVectNorm,numBins,includeScatter,propRegionsRepresented,false);
    title(sprintf('%u superimposed %s gradients: n = %g',numGradients,whatGradients,1/cFree.n));
    % bar(binCenters,propRegionsRepresented)

    % Together:
    % f = figure('color','w'); hold('on')
    % BF_PlotQuantiles(dVectZoom,cgeVect,numBins);
    % BF_PlotQuantiles(dVectZoom,cgeVectNorm,numBins);
    %-------------------------------------------------------------------------------
end


% Plotting:
% f = figure('color','w');

% hold('on');
% for i = 1:length(scalingFactor)
%     plot(xBinCenters,yMeans,'o','color',colors{i})
%     plot(xBinCenters,Exp_3free_fun(xBinCenters,paramStruct{i}),'-','color',colors{i},'LineWidth',2)
% end

%-------------------------------------------------------------------------------
% Plot parametric variation through time:
f = figure('color','w');
subplot(3,1,1);
lambdaEst = arrayfun(@(x)1/paramStruct{x}.n,1:length(paramStruct));
plot(timeVector,lambdaEst,'o-k')
title('Spatial scale')
subplot(3,1,2);
strengthEst = arrayfun(@(x)paramStruct{x}.A,1:length(paramStruct));
plot(timeVector,strengthEst,'o-k')
title('Strength')
subplot(3,1,3);
offSetEst = arrayfun(@(x)paramStruct{x}.B,1:length(paramStruct));
plot(timeVector,offSetEst,'o-k')
title('Offset')
