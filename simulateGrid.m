%-------------------------------------------------------------------------------
% Parameters:
%-------------------------------------------------------------------------------
% Parameters of the 2D spatial coordinate grid:
extentSim = 1/sqrt(2); % physical extent of the grid (square: linear dimension)
resolution = 12; % number of points to split each dimension into
numDims = 3;

% Details of the expression maps to generate (what parameters of what model):
numGradients = 250;
whatGradients = 'spatialLag'; % 'linear' 'poly', 'Gaussian', 'GaussianFixedScale', 'ExpDecaySingle'
ensembleParams = struct();
ensembleParams.rho = 1; % relative strength of spatial autocorrelation (relative to noise)
ensembleParams.d0 = 1; % spatial scale
% Normalization of each gradient:
normalizeHow = 'scaledSigmoid'; % 'zscore','mixedSigmoid','subtractMean'

% Spatial binning for plotting and analysis:
numBins = 25;
includeScatter = false;

% Time Variation:
% --A combination of isotropic spatial scaling, and 'zooming out'--
% spatialScalingFactor = [1,1.5,2,2.5,3]; % in absolute units relative to *original* grid (undistorted by scaling).
spatialScalingFactors = [3.52,5.7,6.6,7.84,10.56,12.8,13.6];
% spatialScalingFactor = [1,1,1,1,1];
% For genes that scale--fraction of average distance
d0scalingFactor = 10;
% This proportion of genes obey the d0 scaling factor, other genes stay put.
propGenesThatScale = 1;
assert(length(d0ScalingFactors)==length(spatialScalingFactors));
numTimePoints = length(spatialScalingFactors);
timeVector = 1:numTimePoints;

%-------------------------------------------------------------------------------
% Now we go through 'time':
f = figure('color','w');
colors = BF_getcmap('spectral',numTimePoints,1);
hold('on')
paramStruct = cell(numTimePoints,1);
for t = timeVector
    % Scale the original map by this much (isotropic stretch):
    scalingFactor = spatialScalingFactors(t);
    % Zoom out from the original map (with the same rules) by this much:
    d0 = d0Values(t);

    %-------------------------------------------------------------------------------
    % Generate the 2-d spatial grid in [X,Y] (on which to work):
    customExtent = GiveMeGridExtent(t);
    [coOrds,X,Y,Z] = MakeGrid(customExtent,resolution,numDims);
    % [coOrds,X,Y,Z] = MakeGrid(extentSim*scalingFactor,resolution,numDims);
    numAreas = length(coOrds);
    % Compute the pairwise distance matrix in the coordinate space:
    dVect = pdist(coOrds,'Euclidean');
    dMat = squareform(dVect);
    fprintf(1,'T = %u, Max distance = %f\n',t,max(dVect));

    %-------------------------------------------------------------------------------
    % Understand the spatial sampling of points across distance bins:
    propRegionsRepresented = AnalyseSpatialSampling(dMat,numAreas,numBins);

    %-------------------------------------------------------------------------------
    % Generate a bunch of spatial gene-expression gradients:
    ensembleParamsT = ensembleParams;
    ensembleParamsT.d0 = mean(dVect)/d0scalingFactor;
    if propGenesThatScale==1
        expData = GenerateEnsemble(whatGradients,numAreas,round(numGradients*propGenesThatScale),dMat,ensembleParamsT);
    else
        fprintf(1,'Only %u%% of genes scale with brain size\n',propGenesThatScale*100);
        % Some genes keep a fixed scale of spatial autocorrelation:
        expDataFixed = GenerateEnsemble(whatGradients,numAreas,round(numGradients*(1 - propGenesThatScale)),dMat,ensembleParams);
        % Genes that scale:
        expDataScale = GenerateEnsemble(whatGradients,numAreas,round(numGradients*propGenesThatScale),dMat,ensembleParamsT);
        expData = [expDataFixed,expDataScale];
    end

    % Normalize each gradient:
    expDataNorm = BF_NormalizeMatrix(expData,normalizeHow);
    % Compute pairwise similarity of gradients as CGE:
    cgeVectNorm = 1 - pdist(expDataNorm,'corr');
    % (unnormalized:)
    cgeVect = 1 - pdist(expData,'corr');

    %-------------------------------------------------------------------------------
    % Plot some individual gradients as color in 2d:
    if t==1 & numDims==2
        f_tmp = figure('color','w');
        set(0, 'currentfigure', f_tmp);  % for figures
        PlotSomeIndividualGradients(expData,coOrds,size(X),false);
        drawnow()
        set(0, 'currentfigure', f);  % for figures
    end

    % %-------------------------------------------------------------------------------
    % % Zoom in on a subset (e.g., for finite-size effects):
    % if zoomFactor > 1
    %     [expDataT,coOrdsT,XT,YT] = ZoomIn(extentSim,zoomFactor,coOrds,expData,X,Y);
    % else
    %     % No zoom at this time point (just copy original variables)
    %     expDataT = expData;
    %     coOrdsT = coOrds;
    %     XT = X;
    %     YT = Y;
    % end

    %-------------------------------------------------------------------------------
    % Rescale to desired brain size:
    % maxExtentDesired = extentSim*scalingFactor;
    % maxExtentCurrent = max(coOrdsT(:));
    % scalingRequired = maxExtentDesired/maxExtentCurrent;
    % % fprintf(1,'T = %u, scaling by %g after zooming by %g\n',t,scalingRequired,zoomFactor);
    % XT = XT*scalingRequired;
    % YT = YT*scalingRequired;
    % coOrdsT = coOrdsT.*scalingRequired;


    %-------------------------------------------------------------------------------
    % Fit:
    [xBinCenters,xThresholds,yMeans,yMedians] = makeQuantiles(dVect,cgeVectNorm,numBins+1);
    [Exp_3free_fun,Stats,paramStruct{t}] = GiveMeFit(xBinCenters',yMeans','exp',true);

    %-------------------------------------------------------------------------------
    % Plot:
    % Unnormalized data:
    % [binCenters,c10,cFree] = PlotWithFit(dVectT,cgeVect,numBins,includeScatter,propRegionsRepresented)
    % title(sprintf('%u superimposed %s gradients: n = %g',numGradients,whatGradients,1/cFree.n))

    % Normalized data:
    subplot(numTimePoints,1,t);
    hold('on')
    [binCenters,c10,cFree] = PlotWithFit(dVect,cgeVectNorm,numBins,includeScatter,propRegionsRepresented,false);
    title(sprintf('d0 %f, scale %f',d0ScalingFactor,scalingFactor));
    % title(sprintf('%u superimposed %s gradients: n = %g',numGradients,whatGradients,1/cFree.n));
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
plot(spatialScalingFactors,lambdaEst,'o-k')
title('Spatial scale')
subplot(3,1,2);
strengthEst = arrayfun(@(x)paramStruct{x}.A,1:length(paramStruct));
plot(spatialScalingFactors,strengthEst,'o-k')
title('Strength')
subplot(3,1,3);
offSetEst = arrayfun(@(x)paramStruct{x}.B,1:length(paramStruct));
plot(spatialScalingFactors,offSetEst,'o-k')
title('Offset')
xlabel('Max Distance')
