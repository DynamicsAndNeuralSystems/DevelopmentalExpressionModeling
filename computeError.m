function sumSquareErrors = computeError(d0,rho)

% Simulate model with specified parameters:
params = GiveMeDefaultParams();
% Over-ride the d0scaling factor and rho values:
params.d0scalingFactor = d0;
params.ensembleParams.rho = rho;
[maxExtent,paramStruct] = simulateGrid(params,false);
lambdaModel = cellfun(@(x)1/x.n,paramStruct);

% Compare to data:
load('parameterFits.mat','paramMeanValues');
lambdaEmpirical = paramMeanValues(1,:)';

sumSquareErrors = sum((lambdaEmpirical-lambdaModel).^2);

end
