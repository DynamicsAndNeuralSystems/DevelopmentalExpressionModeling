function [coOrds,X,Y] = MakeGrid(extent,resolution)
% Generate an (X,Y) grid for embedding spatial maps onto
% --extent is the linear spatial extent in each dimension
% --resolution is the total number of points to split along each dimension (unitless)
%-------------------------------------------------------------------------------

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
