function sumSquareErrors = computeError(d0,rho)

includeF0 = true;

% Simulate model with specified parameters:
params = GiveMeDefaultParams();
% Over-ride the d0scaling factor and rho values:
params.d0scalingFactor = d0;
params.ensembleParams.rho = rho;
[maxExtent,paramStruct] = simulateGrid(params,false);
lambdaModel = cellfun(@(x)1/x.n,paramStruct);
strengthModel = cellfun(@(x)x.A,paramStruct);
f0Model = cellfun(@(x)x.B,paramStruct);

% Compare to data:
load('parameterFits.mat','paramMeanValues');
lambdaEmpirical = paramMeanValues(1,:)';
strengthEmpirical = paramMeanValues(2,:)';
f0Empirical = paramMeanValues(3,:)';

if includeF0
    sumSquareErrors = sum((lambdaEmpirical-lambdaModel).^2) + ...
                sum((strengthEmpirical-strengthModel).^2) + ...
                sum((f0Empirical-f0Model).^2);
else
    sumSquareErrors = sum((lambdaEmpirical-lambdaModel).^2) + ...
                sum((strengthEmpirical-strengthModel).^2);
end

end
