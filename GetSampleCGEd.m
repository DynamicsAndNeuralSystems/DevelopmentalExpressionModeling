
params = GiveMeDefaultParams();
t = 1;

%-------------------------------------------------------------------------------
% Generate the 2d/3d spatial grid in [X,Y] or [X,Y,Z] to work on:
customExtent = GiveMeGridExtent(t);
maxExtent(t) = customExtent(1);
if params.numDims==2
    warning('Using two-dimensions only')
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
% Generate an ensemble of spatial gene-expression gradients:
ensembleParamsT = params.ensembleParams;
ensembleParamsT.d0 = maxExtent(t)/params.d0scalingFactor;
expData = GenerateEnsemble(params.whatGradients,dMat,params.numGradients,ensembleParamsT);

% Normalize each gradient:
expDataNorm = BF_NormalizeMatrix(expData,params.normalizeHow);
% Compute pairwise similarity of gradients as CGE:
cgeVectNorm = 1 - pdist(expDataNorm,'corr');
% (unnormalized:)
% cgeVect = 1 - pdist(expData,'corr');

%-------------------------------------------------------------------------------
% Fit:
[xBinCenters,xThresholds,yMeans,yMedians] = makeQuantiles(dVect,cgeVectNorm,params.numBins+1);
[Exp_3free_fun,Stats,paramStruct{t}] = GiveMeFit(xBinCenters',yMeans','exp',true);

% Normalized data:
f = figure('color','w');
hold('on')
PlotWithFit(dVect,cgeVectNorm,params.numBins,params.includeScatter,...
                                                propRegionsRepresented,false);
title(sprintf('d0 %f',params.d0scalingFactor));
f.Position = [744   826   390   224];
