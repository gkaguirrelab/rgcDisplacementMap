% Validation code
% This script runs a series of validation test on the displacement map in
% order to check if any of the initial principles and assumptions have been
% violated.


%% Set input parameters for displacement map
sampleResolutionDegrees     = 0.01;
maxModeledEccentricity      = 20;
meridianAngleResolutionDeg  = 15;
displacementMapPixelsPerDeg = 7;

%% set range on rand vars that set the target displacement zone end on the meridians
upperLim = 20;
lowerLim = 8;

%% Set number of iteration for bootstrap
numIterations = 8;


%% Make Displacemnt Map

for ii = 1:numIterations
    
    sprintf('Iteration %s start',num2str(ii))
    targetDisplacementAtCardinalMeridiansDeg =(upperLim-lowerLim).*rand(4,1) + lowerLim;
    % verbose must be set to true or else convergenceEccen is not calculated
    [ displacementMapDeg, fitParams, meridianAngles, rgcDisplacementEachMeridian, mRGC_cumulativeEachMeridian, mRF_cumulativeEachMeridian, convergenceEccen(ii,:)] = ...
        makeDisplacementMap('verbose',true,'sampleResolutionDegrees',sampleResolutionDegrees,'maxModeledEccentricity', maxModeledEccentricity, ...
        'meridianAngleResolutionDeg',meridianAngleResolutionDeg,'displacementMapPixelsPerDeg',displacementMapPixelsPerDeg,'verbose',true);
end


%% Mean Disp Zone and SEM

meanDispZone = mean(convergenceEccen);
semDispZone = std(convergenceEccen)./sqrt(numIterations);

figure
hold on

c = categorical({'0','45','90','135'});
c = reordercats(c,{'0','45','90','135'});
bar1 = bar(meanDispZone,'b');
errorbar(meanDispZone,semDispZone,'.r')



