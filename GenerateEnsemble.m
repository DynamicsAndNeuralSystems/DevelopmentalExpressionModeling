function expData = GenerateEnsemble(whatGradients,numAreas,numGradients,dMat,ensembleParams)
% Generates an ensemble of spatial maps, as the matrix expData.
%-------------------------------------------------------------------------------

expData = zeros(numAreas,numGradients);
for i = 1:numGradients
    switch whatGradients
    case 'spatialLag'
        d0 = ensembleParams.d0;
        rho = ensembleParams.rho;
        % d0 = 1; %(max(coOrds(:,1))-min(coOrds(:,1)));
        % rho = 1.5;
        g = GenerateSpatialLagMap(dMat,d0,rho);
    case 'linear'
        X0 = rand*(xRange(2)-xRange(1))+xRange(1);
        Y0 = rand*(yRange(2)-yRange(1))+yRange(1);
        rCoeff = rand(2,1)-0.5;
        g = rCoeff(1)*(X-X0) + rCoeff(2)*(Y-Y0);
    case 'poly'
        pOrderX = randi(1);
        pOrderY = randi(2);
        X0 = rand*(xRange(2)-xRange(1))+xRange(1);
        Y0 = rand*(yRange(2)-yRange(1))+yRange(1);
        rCoeff = rand(2,1)-0.5;
        g = rCoeff(1)*(X-X0).^pOrderX + rCoeff(2)*(Y-Y0).^pOrderY;
    case 'Gaussian'
        A = rand();
        X0 = rand*(xRange(2)-xRange(1))+xRange(1);
        Y0 = rand*(yRange(2)-yRange(1))+yRange(1);
        sigma2X = rand*(xRange(2)-xRange(1))/2;
        sigma2Y = rand*(yRange(2)-yRange(1))/2;
        g = A*exp(-((X-X0).^2./sigma2X + (Y-Y0).^2./sigma2Y));
    case 'GaussianFixedScale'
        A = 1;
        scaleProp = 1/4;
        X0 = rand*(xRange(2)-xRange(1))+xRange(1);
        Y0 = rand*(yRange(2)-yRange(1))+yRange(1);
        sigma2X = (xRange(2)-xRange(1))*scaleProp;
        sigma2Y = (yRange(2)-yRange(1))*scaleProp;
        g = A*exp(-(((X-X0)/sigma2X).^2 + ((Y-Y0)/sigma2Y).^2));
    case 'ExpDecaySingle'
        % d0 = rand*20;
        d0 = ensembleParams.d0;
        rho = ensembleParams.rho;
        % Pick a point as source:
        pointID = randi(length(dMat),1);
        % X0 = rand*(xRange(2)-xRange(1))+xRange(1);
        % Y0 = rand*(yRange(2)-yRange(1))+yRange(1);
        % D = (X-X0).^2 + (Y-Y0).^2;
        distToPoint = dMat(pointID,:);
        g = rho*exp(-distToPoint/d0) + randn(size(distToPoint));
    end
    gStretch = g(:);
    expData(:,i) = gStretch;
end

end
