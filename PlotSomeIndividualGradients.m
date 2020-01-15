function PlotSomeIndividualGradients(expData,coOrds,dims,newFigure)
% Plots some individual gradients as color maps (in unlabeled space)
%-------------------------------------------------------------------------------
if nargin < 4
    newFigure = true;
end
if newFigure
    f = figure('color','w');
end
%-------------------------------------------------------------------------------

ind = randperm(size(expData,2));
for i = 1:12
    ax = subplot(3,4,i);
    f.Position = [1000        1012        1000         326];
    % Plot as a grid:
    imagesc(coOrds(:,1),coOrds(:,2),reshape(expData(:,ind(i)),dims))
    axis('square')
end
colormap(gray)

end
