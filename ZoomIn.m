function [expDataZoom,coOrdsZoom,XZoomed,YZoomed] = ZoomIn(extentSim,zoomFactor,coOrds,expData,X,Y)

% Convert the zoom factor to an extentZoom:
% (such that a zoom factor of 2 would show the middle half of each axis)
extentZoom = extentSim/zoomFactor;

% Get the spatial range to zoom in on:
extentRange = [extentSim/2-extentZoom/2,extentSim/2+extentZoom/2];
isInRange = @(x) x>=extentRange(1) & x<=extentRange(2);
keepMe = isInRange(coOrds(:,1)) & isInRange(coOrds(:,2));

% Take subsets of the relevant variables:
expDataZoom = expData(keepMe,:);
coOrdsZoom = coOrds(keepMe,:);
XKeep = isInRange(X(1,:));
YKeep = isInRange(Y(:,1));
XZoomed = X(YKeep,XKeep);
YZoomed = Y(YKeep,XKeep);
fprintf(1,'We just zoomed in with extent: %u!!!\n',extentZoom);


end
