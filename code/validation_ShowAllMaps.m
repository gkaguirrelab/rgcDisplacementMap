function validation_ShowAllMaps(varargin)

% Transform cone density map to RGC density map


%% Parse input and define variables
p = inputParser;

% Optional anaysis params
p.addParameter('sampleResolutionDegrees',0.01,@isnumeric);
p.addParameter('maxModeledEccentricity',30,@isnumeric);
p.addParameter('targetDisplacementAtCardinalMeridiansDeg',[11 17 17 17],@isnumeric);
p.addParameter('cardinalMeridianAngles',[0 90 180 270],@isnumeric);
p.addParameter('meridianAngleResolutionDeg',15,@isnumeric);
p.addParameter('displacementMapPixelsPerDeg',10,@isnumeric);

% Optional display params
p.addParameter('verbose',false,@islogical);

% parse
p.parse(varargin{:})

close all

%% Setup
% Prepare the regular eccentricity support base
regularSupportPosDeg = ...
    0:p.Results.sampleResolutionDegrees:p.Results.maxModeledEccentricity;

% Get the displacement map
[ ~, fitParams, meridianAngles, rgcDisplacementEachMeridian, mRGC_cumulativeEachMeridian, mRF_cumulativeEachMeridian ] = makeDisplacementMap(  );

%% Loop over the meridians
% Create cone density and RGC density polar maps
for mm = 1:length(meridianAngles)
    
    % obtain cone density
    coneDensityFit = getSplineFitToConeDensity(meridianAngles(mm));
    coneDensitySqDeg = coneDensityFit(regularSupportPosDeg);
    coneDensityEachMeridian(mm,:) = coneDensitySqDeg;

    % obtain the mRF density
    [ mRFDensitySqDeg, mRFtoConeDensityRatio ] = transformConeToMidgetRFDensity( coneDensitySqDeg, 'logitFitParams', fitParams(mm,1:2) );
    mRFDensityEachMeridian(mm,:) = mRFDensitySqDeg;
    mRFtoConeDensityEachMeridian(mm,:) = mRFtoConeDensityRatio;
    
    % obtain the RGC density    
    rgcDensityFit = getSplineFitToRGCDensity(meridianAngles(mm));    
    rgcDensitySqDeg = rgcDensityFit(regularSupportPosDeg);
    rgcDensityEachMeridian(mm,:) = rgcDensitySqDeg;
    
    % obtain the mRGC density
    [ mRGCDensitySqDeg, midgetFraction ] = transformRGCToMidgetRGCDensity( regularSupportPosDeg, rgcDensitySqDeg', 'recipFitParams', fitParams(mm,3:5) );
    mRGCDensityEachMeridian(mm,:) = mRGCDensitySqDeg;
    midgetFractionEachMeridian(mm,:) = midgetFraction;
    
end

imRdim = p.Results.maxModeledEccentricity * p.Results.displacementMapPixelsPerDeg * 2;

% Show the maps
figure; imagesc(convertPolarMapToImageMap(coneDensityEachMeridian, imRdim));
figure; imagesc(convertPolarMapToImageMap(mRFDensityEachMeridian, imRdim));
figure; imagesc(convertPolarMapToImageMap(mRFtoConeDensityEachMeridian, imRdim));
figure; imagesc(convertPolarMapToImageMap(rgcDensityEachMeridian, imRdim));
figure; imagesc(convertPolarMapToImageMap(mRGCDensityEachMeridian, imRdim));
figure; imagesc(convertPolarMapToImageMap(midgetFractionEachMeridian, imRdim));
figure; imagesc(convertPolarMapToImageMap(rgcDisplacementEachMeridian, imRdim));
figure; imagesc(convertPolarMapToImageMap(mRGC_cumulativeEachMeridian, imRdim));
figure; imagesc(convertPolarMapToImageMap(mRF_cumulativeEachMeridian, imRdim));

% create a

end % function


function imageMap = convertPolarMapToImageMap(polarMap, imRdim)
maxMapDensity = max(polarMap(:));
imP=polarMap'./maxMapDensity;
imR = PolarToIm (imP, 0, 1, imRdim, imRdim);
imageMap = imrotate(imR .* maxMapDensity,-90);
end

% function displacedValues = applyDisplacement(originalValues, displacementInDeg, regularSupportPosDeg)
% 
% countsPerRing = calcCumulative(regularSupportPosDeg, originalValues);
% displacedValues=zeros(size(countsPerRing));
% [~, displacedIdx] = arrayfun(@(x) min(abs(x-regularSupportPosDeg)), regularSupportPosDeg-displacementInDeg);
% displacedIdx = displacedIdx-1;
% targetList = unique(displacedIdx);
% targetList = targetList(targetList > 0);
% piled = arrayfun(@(x) sum(countsPerRing(displacedIdx==targetList(x))), 1:1:length(targetList));
% displacedValues(targetList)=piled;
% end
