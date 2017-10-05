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

imRdim = (length(xSmps) * 2)-1;

maxMapDensity = max(rgcDensityMatrix(:));
imP=rgcDensityMatrix'./maxMapDensity;
imR = PolarToIm (imP, 0, 1, imRdim, imRdim);
rgcDensityMapDeg = imrotate(imR .* maxMapDensity,-90);





%% Get displacement map
[ ~, fitParams, meridianAngles, rgcDisplacementEachMeridian, mRGC_cumulativeEachMeridian, mRF_cumulativeEachMeridian ] = ...
    makeDisplacementMap(...
    'maxModeledEccentricity', maxModeledEccentricity, ...
    'displacementMapPixelsPerDeg', displacementMapPixelsPerDeg, ...
    'verbose',true);


maxMapDensity = max(rgcDisplacementEachMeridian(:));
imP=rgcDisplacementEachMeridian'./maxMapDensity;
imR = PolarToIm (imP, 0, 1, imRdim, imRdim);
displacementMapDeg = imrotate(imR .* maxMapDensity,-90);


%% Apply warp

smps = linspace(-1*maxModeledEccentricity,maxModeledEccentricity,(length(xSmps)*2)-1);

[sampleBaseX,sampleBaseY] = meshgrid(smps,smps);
warpedRGC = applyDisplacementMap( rgcDensityMapDeg, displacementMapDeg, sampleBaseX, sampleBaseY );

%% fix tears in warped image
%// identify indices valid for the warpedRGC
idxgood=~(isnan(sampleBaseX) | isnan(sampleBaseY) | isnan(warpedRGC)); 

%// re-interpolate scattered data over the sampleBase grid
warpedRGC_filled = griddata( sampleBaseX(idxgood),sampleBaseY(idxgood),warpedRGC(idxgood), sampleBaseX, sampleBaseY ) ;


%% Smooth the final map
H = fspecial('gaussian',7,2);
warpedRGC_gaussian = imfilter(warpedRGC_filled,H,'replicate'); 

figure;
surf(warpedRGC_gaussian)






