function PlotSomeIndividualGradients(expData,dims)

f = figure('color','w');
for i = 1:12
    subplot(3,4,i)
    f.Position = [1000        1012        1000         326];
    % Plot as a grid:
    imagesc(reshape(expDataZoom(:,i),size(XZoomed)))
    axis('square')
end
colormap(gray)

end
