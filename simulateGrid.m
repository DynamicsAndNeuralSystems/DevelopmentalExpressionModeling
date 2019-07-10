% Make a grid in 2d space
whatGradients = 'GaussianFixedScale'; % 'linear' 'poly', 'Gaussian', 'GaussianFixedScale'
numGradients = 40;
xRange = [0,0.2,0.6,1,1.1,1.4,2.5,3.1,4,4.3,6,6.9,8.1,9.3,10];
yRange = linspace(0,10,15);

% Distance (m) corresponding to a unit change in each axis
dScale = 1;

[X,Y] = meshgrid(xRange,yRange);
coOrds = [X(:),Y(:)];
numAreas = length(coOrds);

% Compute gradients through the space:
expData = zeros(numAreas,numGradients);
for i = 1:numGradients
    switch whatGradients
    case 'linear'
        X0 = rand*(xRange(end)-xRange(1))+xRange(1);
        Y0 = rand*(yRange(end)-yRange(1))+yRange(1);
        rCoeff = rand(2,1)-0.5;
        g = rCoeff(1)*(X-X0) + rCoeff(2)*(Y-Y0);
    case 'poly'
        pOrderX = randi(1);
        pOrderY = randi(2);
        X0 = rand*(xRange(end)-xRange(1))+xRange(1);
        Y0 = rand*(yRange(end)-yRange(1))+yRange(1);
        rCoeff = rand(2,1)-0.5;
        g = rCoeff(1)*(X-X0).^pOrderX + rCoeff(2)*(Y-Y0).^pOrderY;
    case 'Gaussian'
        A = rand();
        X0 = rand*(xRange(end)-xRange(1))+xRange(1);
        Y0 = rand*(yRange(end)-yRange(1))+yRange(1);
        sigma2X = rand*(xRange(end)-xRange(1))/2;
        sigma2Y = rand*(yRange(end)-yRange(1))/2;
        g = A*exp(-((X-X0).^2./sigma2X + (Y-Y0).^2./sigma2Y));
    case 'GaussianFixedScale'
        A = 1;
        X0 = rand*(xRange(end)-xRange(1))+xRange(1);
        Y0 = rand*(yRange(end)-yRange(1))+yRange(1);
        sigma2X = (xRange(end)-xRange(1))/5;
        sigma2Y = (yRange(end)-yRange(1))/5;
        g = A*exp(-((X-X0).^2./sigma2X + (Y-Y0).^2./sigma2Y));
    end
    gStretch = g(:);
    expData(:,i) = gStretch;
end
expDataZ = zscore(expData);

% Compute pairwise distance matrix:
dVect = pdist(dScale*coOrds,'Euclidean');

% Compute pairwise similarity of gradient:
cgeVect = 1 - pdist(expDataZ,'corr');

%-------------------------------------------------------------------------------
% Plot against each other:
f = figure('color','w');
f.Position = [1000        1012        1000         326];
% Plot as a grid:
subplot(1,2,1)
surf(X,Y,reshape(sum(expDataZ,2),size(X)))
xlabel('x')
ylabel('y')
zlabel('z')
%-------------------------------------------------------------------------------
subplot(1,2,2)
plot(dVect,cgeVect,'.k')
xlabel('Distance (m)')
ylabel('CGE')
title(sprintf('%u superimposed %s gradients',numGradients,whatGradients))
