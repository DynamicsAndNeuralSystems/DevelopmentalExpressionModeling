function [maxExtent,paramStruct] = simulateGrid(params,doPlot)
%-------------------------------------------------------------------------------
if nargin < 1
    params = GiveMeDefaultParams();
end
if nargin < 2
    doPlot = false;
end
%-------------------------------------------------------------------------------
numTimePoints = 7; %length(spatialScalingFactors);
timeVector = 1:numTimePoints;

%-------------------------------------------------------------------------------
% Now we go through 'time':
if doPlot
    f = figure('color','w');
    colors = BF_getcmap('spectral',numTimePoints,1);
    hold('on')
end
paramStruct = cell(numTimePoints,1);
maxExtent = zeros(numTimePoints,1);
for t = timeVector
    % Scale the original map by this much (isotropic stretch):
    % scalingFactor = spatialScalingFactors(t);
    % Zoom out from the original map (with the same rules) by this much:
    % d0 = d0Values(t);

    %-------------------------------------------------------------------------------
    % Generate the 2-d spatial grid in [X,Y] (on which to work):
    customExtent = GiveMeGridExtent(t);
    maxExtent(t) = customExtent(1);
    if params.numDims==2
        customExtent = customExtent(1:2);
    end
    [coOrds,X,Y,Z] = MakeGrid(customExtent,params.resolution,params.numDims,params.subsampleSpace);
    % [coOrds,X,Y,Z] = MakeGrid(extentSim*scalingFactor,resolution,numDims);
    numAreas = length(coOrds);
    % Compute the pairwise distance matrix in the coordinate space:
    dVect = pdist(coOrds,'Euclidean');
    dMat = squareform(dVect);
    fprintf(1,'T = %u, Max distance = %f\n',t,max(dVect));

    %-------------------------------------------------------------------------------
    % Understand the spatial sampling of points across distance bins:
    propRegionsRepresented = AnalyseSpatialSampling(dMat,numAreas,params.numBins);

    %-------------------------------------------------------------------------------
    % Generate a bunch of spatial gene-expression gradients:
    ensembleParamsT = params.ensembleParams;
    ensembleParamsT.d0 = mean(dVect)/params.d0scalingFactor;
    if params.propGenesThatScale==1
        expData = GenerateEnsemble(params.whatGradients,dMat,round(params.numGradients*params.propGenesThatScale),ensembleParamsT);
    else
        fprintf(1,'Only %u%% of genes scale with brain size\n',params.propGenesThatScale*100);
        % Some genes keep a fixed scale of spatial autocorrelation:
        expDataFixed = GenerateEnsemble(params.whatGradients,dMat,round(params.numGradients*(1 - params.propGenesThatScale)),ensembleParams);
        % Genes that scale:
        expDataScale = GenerateEnsemble(params.whatGradients,dMat,round(params.numGradients*params.propGenesThatScale),ensembleParamsT);
        expData = [expDataFixed,expDataScale];
    end

    % Normalize each gradient:
    expDataNorm = BF_NormalizeMatrix(expData,params.normalizeHow);
    % Compute pairwise similarity of gradients as CGE:
    cgeVectNorm = 1 - pdist(expDataNorm,'corr');
    % (unnormalized:)
    % cgeVect = 1 - pdist(expData,'corr');

    %-------------------------------------------------------------------------------
    % Plot some individual gradients as color in 2d:
    if t==1 & params.numDims==2 & isempty(params.subsampleSpace)
        f_tmp = figure('color','w');
        set(0, 'currentfigure', f_tmp);  % for figures
        PlotSomeIndividualGradients(expData,coOrds,size(X),true);
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
    [xBinCenters,xThresholds,yMeans,yMedians] = makeQuantiles(dVect,cgeVectNorm,params.numBins+1);
    [Exp_3free_fun,Stats,paramStruct{t}] = GiveMeFit(xBinCenters',yMeans','exp',true);

    %-------------------------------------------------------------------------------
    % Plot:
    % Unnormalized data:
    % [binCenters,c10,cFree] = PlotWithFit(dVectT,cgeVect,numBins,includeScatter,propRegionsRepresented)
    % title(sprintf('%u superimposed %s gradients: n = %g',numGradients,params.whatGradients,1/cFree.n))

    % Normalized data:
    if doPlot
        subplot(numTimePoints,1,t);
        hold('on')
        [binCenters,c10,cFree] = PlotWithFit(dVect,cgeVectNorm,params.numBins,params.includeScatter,...
                                            params.propRegionsRepresented,false);
        title(sprintf('d0 %f',params.d0ScalingFactor));
    end
    % title(sprintf('%u superimposed %s gradients: n = %g',numGradients,params.whatGradients,1/cFree.n));
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
lambdaEst = cellfun(@(x)1/x.n,paramStruct);
strengthEst = cellfun(@(x)x.A,paramStruct);
offSetEst = cellfun(@(x)x.B,paramStruct);

if doPlot
    f = figure('color','w');
    subplot(3,1,1);
    plot(maxExtent,lambdaEst,'o-k')
    title('Spatial scale')
    subplot(3,1,2);
    plot(maxExtent,strengthEst,'o-k')
    title('Strength')
    subplot(3,1,3);
    plot(maxExtent,offSetEst,'o-k')
    title('Offset')
    xlabel('Max Distance')
end
end
