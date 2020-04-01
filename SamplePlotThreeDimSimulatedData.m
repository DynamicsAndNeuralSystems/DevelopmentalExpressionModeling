
% The time point:
t = 1;

params = GiveMeDefaultParams();
params.numGradients = 5;
params.ensembleParams.rho = 0.3;
params.resolution = 25; % number of points to split each dimension into
params.subsampleSpace = [];
customExtent = GiveMeGridExtent(t);
maxExtent(t) = customExtent(1);
[coOrds,X,Y,Z] = MakeGrid(customExtent,params.resolution,params.numDims,params.subsampleSpace);
dVect = pdist(coOrds,'Euclidean');
dMat = squareform(dVect);
ensembleParamsT = params.ensembleParams;
ensembleParamsT.d0 = maxExtent(t)/params.d0scalingFactor;
expData = GenerateEnsemble(params.whatGradients,dMat,params.numGradients,ensembleParamsT);


f = figure('color','w');
whatMap = 2;
theMap = reshape(expData(:,whatMap),[params.resolution,params.resolution,params.resolution]);
for i = 1:params.resolution
    subplot(5,5,i)
    theMapSlice = squeeze(theMap(i,:,:));
    imagesc(1:params.resolution,1:params.resolution,theMapSlice)
    axis('square')
end
colormap(BF_getcmap('redyellowblue',5,false,true))
