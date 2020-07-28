function sumSquareErrors = computeError(d0,rho)
% Get error of model-data fit to scaling of CGE(d) parameters
%-------------------------------------------------------------------------------

includeF0 = false;

% Simulate model with specified parameters:
params = GiveMeDefaultParams();
% Over-ride the d0scaling factor and rho values:
params.d0scalingFactor = d0;
params.ensembleParams.rho = rho;
[maxExtent,paramStruct] = simulateGrid(params,false);
lambdaModel = cellfun(@(x)1/x.n,paramStruct);
strengthModel = cellfun(@(x)x.A,paramStruct);
f0Model = cellfun(@(x)x.B,paramStruct);

%-------------------------------------------------------------------------------
% Load results from fitted data:
%-------------------------------------------------------------------------------
params = struct();
params.doSubsample = true;
params.thisBrainDiv = 'brain';
params.thisCellType = 'allCellTypes';
params.includeAdult = false;
[params,fittedParams,CIs,goodTimePoint] = LoadParameterFits(params);

lambdaEmpirical = cellfun(@(x)1/x.n,fittedParams);
strengthEmpirical = cellfun(@(x)x.A,fittedParams);
f0Empirical = cellfun(@(x)x.B,fittedParams);

if includeF0
    sumSquareErrors = sum(abs((lambdaEmpirical-lambdaModel)./lambdaEmpirical)) + ...
                sum(abs((strengthEmpirical-strengthModel)./strengthEmpirical)) + ...
                sum(abs((f0Empirical-f0Model)./f0Empirical));
else
    sumSquareErrors = sum(abs((lambdaEmpirical-lambdaModel)./lambdaEmpirical)) + ...
                sum(abs((strengthEmpirical-strengthModel)./strengthEmpirical));
end

end
