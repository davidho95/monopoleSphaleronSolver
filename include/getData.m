

% Script to scrape and postprocess output data from C++ monopole sphaleron 
% gradient descent solver

bigNum = 1e10;
fileid = fopen('coords.txt', 'r');
coordsList = fscanf(fileid, '%f %f %f', [3 bigNum])';
[coordsList, sortIdx] = sortrows(coordsList, 3);
matSize = max(coordsList) + 1;
xPoints = reshape(coordsList(:,1), matSize);
yPoints = reshape(coordsList(:,2), matSize);
zPoints = reshape(coordsList(:,3), matSize);


fileid = fopen('higgsData.txt');
higgsData = fscanf(fileid, '%f %f %f', [1 prod(matSize)])';
higgsData = higgsData(sortIdx, :);
higgsCube = reshape(higgsData', [matSize]);
% higgsNorms = squeeze(vecnorm(higgsCube));

fileid = fopen('magneticFieldData.txt');
magneticField = fscanf(fileid, '%f %f %f', [3 prod(matSize)])';
magneticField = magneticField(sortIdx, :);
magX = reshape(magneticField(:,1), matSize);
magY = reshape(magneticField(:,2), matSize);
magZ = reshape(magneticField(:,3), matSize);

% fileid = fopen('gaugeData.txt');
% gaugeField = fscanf(fileid, '%f %f %f', [3 prod(matSize)])';
% gaugeField = gaugeField(sortIdx, :);
% gaugeX = reshape(gaugeField(:,1), matSize);
% gaugeY = reshape(gaugeField(:,2), matSize);
% gaugeZ = reshape(gaugeField(:,3), matSize);

fileid = fopen('energyData.txt');
energyDensity = fscanf(fileid, '%f %f %f', prod(matSize))';
energyDensity = energyDensity(sortIdx);
energyDensity = reshape(energyDensity, matSize);

% fileid = fopen('gradData.txt');
% gradSq = fscanf(fileid, '%f %f %f', prod(matSize))';
% gradSq = gradSq(sortIdx);
% gradSq = reshape(gradSq, matSize);

divB = circshift(magX, -1, 1) + circshift(magY, -1, 2) + circshift(magZ, -1, 3) - magX - magY - magZ;