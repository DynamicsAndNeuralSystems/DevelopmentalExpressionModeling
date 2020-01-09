function [expDataZoom,coOrdsZoom,XZoomed,YZoomed] = ZoomIn(extentSim,extentZoom,coOrds,expData,X,Y)

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
