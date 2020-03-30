%-------------------------------------------------------------------------------
% Plots for schematic:
%-------------------------------------------------------------------------------

% Generate ensemble:
[coOrds,X,Y,Z] = MakeGrid([20,20],50,2,[]);
numAreas = length(coOrds);
dVect = pdist(coOrds,'Euclidean');
dMat = squareform(dVect);
ensembleParams = struct();
ensembleParams.rho = 1; % relative strength of spatial autocorrelation (relative to noise)
ensembleParams.d0 = 4; % spatial scale
expData = GenerateEnsemble('spatialLag',dMat,20,ensembleParams);

%-------------------------------------------------------------------------------
% Understand the spatial sampling of points across distance bins:
PlotSomeIndividualGradients(expData,coOrds,size(X),true)
