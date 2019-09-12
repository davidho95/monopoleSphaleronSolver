function fig = magneticFieldQuiver(vev, g, B, L, baseDir, varargin)
%Produce 

dirString = strcat(baseDir, '/saddleData', strrep(num2str(vev), '.', '_'));
cd(dirString)
getData;

magX0 = magX + 4*pi*B / (g * L^2);
magNorms = sqrt(magX0.^2 + magY.^2 + magZ.^2);
magXDir = magX0 ./ magNorms;
magYDir = magY ./ magNorms;
magZDir = magZ ./ magNorms;

xPointsScaled = sqrt(2)*g*vev*(xPoints - matSize(1)/2 + 0.5);
yPointsScaled = sqrt(2)*g*vev*(yPoints - matSize(2)/2 + 0.5);
zPointsScaled = sqrt(2)*g*vev*(zPoints - matSize(3)/2 + 0.5);

delta = 1;
if numel(varargin) >= 1
    delta = varargin{1};
end

newFig = true;
if numel(varargin) >= 2
    newFig = varargin{2};
end

if newFig
    fig = figure();
else
    fig = gcf;
end

zSlicePoint = matSize(3) / 2 + 2;
magYDir = magYDir(1:delta:end, 1:delta:end, zSlicePoint);
magXDir = magXDir(1:delta:end, 1:delta:end, zSlicePoint);

xSlice = xPointsScaled(1:delta:end,1:delta:end,zSlicePoint);
ySlice = yPointsScaled(1:delta:end,1:delta:end,zSlicePoint);

quiver(ySlice, xSlice, magYDir, magXDir, 'k')

fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

axis('tight');

end

