function fig = energyDensityPlot(vev, g, B, L, baseDir, varargin)
%Produce energy density contour plot

dirString = strcat(baseDir, '/saddleData', strrep(num2str(vev), '.', '_'));
cd(dirString)
getData;

eDensity0 = energyDensity - 8*pi^2*B^2 / (g^2*L^4);
eDensityScaled = eDensity0 / (g*vev)^4;

xPointsScaled = sqrt(2)*g*vev*(xPoints - matSize(1)/2 + 0.5);
yPointsScaled = sqrt(2)*g*vev*(yPoints - matSize(2)/2 + 0.5);
zPointsScaled = sqrt(2)*g*vev*(zPoints - matSize(3)/2 + 0.5);

fig = figure;

delta = 1;
if numel(varargin) >= 1
    delta = varargin{1};
end

xSlice = xPointsScaled(1:delta:end,1,1);
ySlice = yPointsScaled(1,1:delta:end,1);
zSlice = squeeze(zPointsScaled(1,1,1:delta:end));

minFrac = 0;
if numel(varargin) >= 2
    minFrac = varargin{2};
end

numLvls = 10;
if numel(varargin) >= 3
    numLvls = varargin{3};
end

newFig = true;
if numel(varargin) >= 4
    newFig = varargin{4};
end

if newFig
    fig = figure();
else
    fig = gcf;
end
    
eMax = max(eDensityScaled(:));
eMin = eMax*minFrac;
lvls = linspace(eMin, eMax, numLvls);
contourslice(yPointsScaled, xPointsScaled, zPointsScaled, eDensityScaled,...
    xSlice, ySlice, zSlice, lvls);

% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];

end

