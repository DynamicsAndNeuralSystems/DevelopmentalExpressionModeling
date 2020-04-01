function params = GiveMeDefaultParams()

%-------------------------------------------------------------------------------
% Parameters of the 2D/3D spatial coordinate grid:
%-------------------------------------------------------------------------------
params.resolution = 50; % number of points to split each dimension into
params.numDims = 3;
params.extentSim = 1/sqrt(params.numDims); % physical extent of the grid (square: linear dimension)

%-------------------------------------------------------------------------------
% Details of the expression maps to generate (what parameters of what generative model):
%-------------------------------------------------------------------------------
params.numGradients = 1861;
params.whatGradients = 'spatialLag'; % 'linear' 'poly', 'Gaussian', 'GaussianFixedScale', 'ExpDecaySingle'
params.ensembleParams = struct();
params.ensembleParams.rho = 0.16; % relative strength of spatial autocorrelation (relative to noise)
params.ensembleParams.d0 = []; % spatial scale

% Normalization of each gradient:
params.normalizeHow = 'scaledSigmoid'; % 'zscore','mixedSigmoid','subtractMean'

%-------------------------------------------------------------------------------
% Spatial binning for plotting and analysis:
%-------------------------------------------------------------------------------
params.subsampleSpace = 1000; % only keep this many random points for spatial analysis
params.numBins = 20;
params.includeScatter = false;

%-------------------------------------------------------------------------------
% Time Variation:
% --A combination of isotropic spatial scaling, and 'zooming out'--
% spatialScalingFactor = [1,1.5,2,2.5,3]; % in absolute units relative to *original* grid (undistorted by scaling).
% spatialScalingFactors = [3.52,5.7,6.6,7.84,10.56,12.8,13.6];
% spatialScalingFactor = [1,1,1,1,1];
% For genes that scale--fraction of average distance
params.d0scalingFactor = 8.41;
% This proportion of genes obey the d0 scaling factor, other genes stay put.
params.propGenesThatScale = 1;

end
