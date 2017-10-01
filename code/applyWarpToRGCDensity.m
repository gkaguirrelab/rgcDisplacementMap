%% applyWarpToRGC
% 
% This script generates a 2D RGC map from looping over meridian angles with
% getSplineFitToRGCDensity.m and applies the displacement map from 
% makeDisplacementMap.m withte function warpImage.m
%


%% Set up params

maxModeledEccentricity      = 20;
displacementMapPixelsPerDeg = 10;
mapResolutionDegrees        = 1/displacementMapPixelsPerDeg;
xSmps = mapResolutionDegrees:mapResolutionDegrees:maxModeledEccentricity;


%% make RGC Density Map

angularResolutionDegrees = 15;
angles = 0:angularResolutionDegrees:360-angularResolutionDegrees;

for ii = 1:length(angles)
    sprintf('Fitting %s degree meridian: %s% Done',num2str(ii),num2str(round(ii/length(angles),2)))
    rgcDensityFit          = getSplineFitToRGCDensity(angles(ii));
    rgcDensityMatrix(ii,:) = rgcDensityFit(xSmps);
end

imRdim = length(xSmps) * 2;
maxDensityDeg = max(rgcDensityMatrix(:));
imP=rgcDensityMatrix'./maxDensityDeg;
imR = PolarToIm (imP, 0, 1, imRdim, imRdim);
rgcDensityMapDeg = imrotate(imR .* maxDensityDeg,-90);


%% Get displacement map

[ displacementMapDeg, fitParams, meridianAngles, rgcDisplacementEachMeridian, mRGC_cumulativeEachMeridian, mRF_cumulativeEachMeridian, convergenceEccen] = ...
    makeDisplacementMap('verbose',true,'maxModeledEccentricity',maxModeledEccentricity, 'displacementMapPixelsPerDeg',displacementMapPixelsPerDeg);


%% Apply warp

smps = linspace(-1*maxModeledEccentricity,maxModeledEccentricity,length(xSmps)*2);

[sampleBaseX,sampleBaseY] = meshgrid(smps,smps);
warpedRGC = warpImage( rgcDensityMapDeg, displacementMapDeg, sampleBaseX, sampleBaseY );