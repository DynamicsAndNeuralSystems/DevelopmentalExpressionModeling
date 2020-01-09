function [coOrds,X,Y] = MakeGrid(extent,resolution)

if nargin < 1
    extent = 10;
end
if nargin < 2
    resolution = 25;
end

%-------------------------------------------------------------------------------
xRange = linspace(0,extent,resolution);
yRange = linspace(0,extent,resolution);

[X,Y] = meshgrid(xRange,yRange);
coOrds = [X(:),Y(:)];

end
