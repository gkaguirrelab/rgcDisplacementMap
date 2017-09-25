% Validation code 



%% Set input parameters for displacement map
sampleResolutionDegrees     = 0.01;
maxModeledEccentricity      = 30;
meridianAngleResolutionDeg  = 15;
displacementMapPixelsPerDeg = 10;

%% Make Displacemnt Map
[ displacementMapDeg, fitParams, meridianAngles, rgcDisplacementEachMeridian, mRGC_cumulativeEachMeridian, mRF_cumulativeEachMeridian ] = ...
    makeDisplacementMap('verbose',true,'sampleResolutionDegrees',sampleResolutionDegrees,'maxModeledEccentricity', maxModeledEccentricity, ... 
    'meridianAngleResolutionDeg',meridianAngleResolutionDeg,'displacementMapPixelsPerDeg',displacementMapPixelsPerDeg);

