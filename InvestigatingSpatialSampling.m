%-------------------------------------------------------------------------------
% Visualize points at each distance bin:
numBins = 20;
doPlot = false;
anglesOrSpace = 'space';

dMatUpper = dMat;
dMatUpper(tril(true(size(dMat)))) = 0;
xThresholds = arrayfun(@(x)quantile(dMatUpper(dMatUpper>0),x),linspace(0,1,numBins+1));
xThresholds(end) = xThresholds(end) + eps; % make sure all data included in final bin
f = figure('color','w');
propRegionsRepresented = zeros(numBins,1);
for i = 1:numBins
    subplot(5,4,i)
    axis('square')

    % Get points in bin index i:
    [isInBin_i,isInBin_j] = find(dMatUpper>xThresholds(i) & dMatUpper<=xThresholds(i+1));
    regionIsRepresented = union(isInBin_i,isInBin_j);
    propRegionsRepresented(i) = length(regionIsRepresented)/numAreas;

    switch anglesOrSpace
    case 'space'
        % Plot points contributing to the distance bin in coordinate space:
        plot(coOrds(regionIsRepresented,1),coOrds(regionIsRepresented,2),'.r')

    case 'angles'
        % Get distribution of angles:
        numPointsInBin = length(isInBin_i);
        angles = zeros(numPointsInBin,1);
        for j = 1:numPointsInBin
            % Pointer from 1 -> 2:
            diffVector = coOrds(isInBin_i(j),:) - coOrds(isInBin_j(j),:);
            % Determine angle from inverse tan:
            angles(j) = atan(diffVector(2)/diffVector(1));
        end
        polarhistogram(abs(angles))
    end
end
