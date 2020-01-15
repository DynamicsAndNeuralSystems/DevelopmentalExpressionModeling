function [coOrds,X,Y,Z] = MakeGrid(extent,resolution,numDims,subsampleSpace)
% Generate an (X,Y,[Z]) grid for embedding spatial maps onto
% --extent is the linear spatial extent in each dimension
% --resolution is the total number of points to split along each dimension (unitless)
%-------------------------------------------------------------------------------

if nargin < 1
    extent = 10;
end
if nargin < 2
    resolution = 25;
end
if nargin < 3
    numDims = 2;
end
assert(ismember(numDims,[2,3]));
if length(extent)==1
    % Isotropic by defualt
    extent = ones(numDims,1)*extent;
end
assert(numDims==length(extent))
if nargin < 4
    subsampleSpace = []; % keep all samples
end

%-------------------------------------------------------------------------------
% Set up isotropic numDims-dimensional grid:
xRange = linspace(0,extent(1),resolution);
yRange = linspace(0,extent(2),resolution);
if numDims==3
    zRange = linspace(0,extent(3),resolution);
end

%-------------------------------------------------------------------------------
if numDims==2
    [X,Y] = meshgrid(xRange,yRange);
    coOrds = [X(:),Y(:)];
    Z = [];
else
    [X,Y,Z] = meshgrid(xRange,yRange,zRange);
    coOrds = [X(:),Y(:),Z(:)];
end

%-------------------------------------------------------------------------------
% Sub-sample
if ~isempty(subsampleSpace)
    rp = randperm(length(coOrds));
    keepMe = rp(1:subsampleSpace);
    coOrds = coOrds(keepMe,:);
    fprintf(1,'Subsampling %u down to %u spatial points\n',length(rp),length(coOrds));
end

end
