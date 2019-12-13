function [coOrds,X,Y] = MakeGrid(extent,resolution)

if nargin < 1
    extent = 10;
end
if nargin < 2
    resolution = 25;
end

% xRange = [0,0.2,0.6,1,1.1,1.4,2.5,3.1,4,4.3,6,6.9,8.1,9.3,10];
xRange = linspace(0,extent,resolution);
yRange = linspace(0,extent,resolution);

[X,Y] = meshgrid(xRange,yRange);
coOrds = [X(:),Y(:)];

end
