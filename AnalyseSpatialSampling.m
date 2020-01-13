function propRegionsRepresented = AnalyseSpatialSampling(dMat,numAreas,numBins)
% Spatial sampling by distance bin is uneven
% This function makes this clear.
% Takes in the distance matrix and the number of distance bins used to partition.
%-------------------------------------------------------------------------------

dVect = squareform(dMat);
dMatUpper = dMat;
dMatUpper(tril(true(size(dMat)))) = 0;

% Threshold distances:
[~,xThresholds] = makeQuantiles(dVect,[],numBins+1);

% Compute the proportion of coordinates included in each bin:
propRegionsRepresented = zeros(numBins,1);
for i = 1:numBins
    % Get points in bin index i:
    [isInBin_i,isInBin_j] = find(dMatUpper > xThresholds(i) & dMatUpper <= xThresholds(i+1));
    regionIsRepresented = union(isInBin_i,isInBin_j);
    propRegionsRepresented(i) = length(regionIsRepresented)/numAreas;
end

end
