function fig = higgsSurface(vev, g, baseDir)
%Produce higgs surface plot

dirString = strcat(baseDir, '/saddleData', strrep(num2str(vev), '.', '_'));
cd(dirString)
getData;

xPointsScaled = sqrt(2)*g*vev*(xPoints - matSize(1)/2 + 0.5);
yPointsScaled = sqrt(2)*g*vev*(yPoints - matSize(2)/2 + 0.5);

fig = figure;

zSlicePoint = matSize(3) / 2 + 2;

xSlice = xPointsScaled(:,:,zSlicePoint);
ySlice = yPointsScaled(:,:,zSlicePoint);

higgsCubeScaled = higgsCube / (vev / sqrt(2));

surf(ySlice, xSlice, higgsCubeScaled(:,:,zSlicePoint))

fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

zlim([0 1])

axis('tight');

end

