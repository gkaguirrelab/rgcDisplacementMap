% Validation code 
% This script runs a series of validation test on the displacement map in
% order to check if any of the initial principles and assumptions have been
% violated.
%
% First Check:  


%% Set input parameters for displacement map
sampleResolutionDegrees     = 0.01;
maxModeledEccentricity      = 30;
meridianAngleResolutionDeg  = 15;
displacementMapPixelsPerDeg = 10;

%% Make Displacemnt Map
[ displacementMapDeg, fitParams, meridianAngles, rgcDisplacementEachMeridian, mRGC_cumulativeEachMeridian, mRF_cumulativeEachMeridian ] = ...
    makeDisplacementMap('verbose',true,'sampleResolutionDegrees',sampleResolutionDegrees,'maxModeledEccentricity', maxModeledEccentricity, ... 
    'meridianAngleResolutionDeg',meridianAngleResolutionDeg,'displacementMapPixelsPerDeg',displacementMapPixelsPerDeg);


[rgcDiffDisplacementEachMeridian,diffDisplacementMapDeg] = calcAndPlotDiffDisplacement(rgcDisplacementEachMeridian,'maxModeledEccentricity', ...
    maxModeledEccentricity, 'displacementMapPixelsPerDeg',displacementMapPixelsPerDeg, 'verbose', true);


