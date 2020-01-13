function PlotSomeIndividualGradients(expData,coOrds,dims)
% Plots some individual gradients as color maps (in unlabeled space)
%-------------------------------------------------------------------------------

f = figure('color','w');
for i = 1:12
    ax = subplot(3,4,i);
    f.Position = [1000        1012        1000         326];
    % Plot as a grid:
    imagesc(coOrds(:,1),coOrds(:,2),reshape(expData(:,i),dims))
    axis('square')
end
colormap(gray)

end
