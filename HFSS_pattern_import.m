function [pattern_phitheta, phi, theta] = HFSS_pattern_import(relative_filename)

patternData = csvread(relative_filename); % import csv

% Extract phi/theta values from custom pattern

chktheta = (patternData(:,2) == patternData(1,2));
blockLen = length(chktheta(chktheta~=0));
nCols = size(patternData,1)/blockLen;
thetaBlocks = reshape(patternData(:,2),blockLen,nCols);
phiBlocks = reshape(patternData(:,1),blockLen,nCols);

theta = thetaBlocks(1,:);
phi = phiBlocks(:,1).';

pattern_phitheta = reshape(patternData(:,3),blockLen,nCols).';
