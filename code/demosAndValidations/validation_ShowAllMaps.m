function validation_ShowAllMaps(varargin)

% Transform cone density map to RGC density map


%% Parse input and define variables
p = inputParser; p.KeepUnmatched = true;

% Optional anaysis params
p.addParameter('sampleResolutionDegrees',0.01,@isnumeric);
p.addParameter('maxModeledEccentricity',30,@isnumeric);
p.addParameter('meridianAngleResolutionDeg',1,@isnumeric);
p.addParameter('displacementMapPixelsPerDeg',10,@isnumeric);
p.addParameter('subjectName', 'reportedAverage', @ischar);

% Optional display and ouput params
p.addParameter('verbose',true,@islogical);
p.addParameter('savePlots',true,@islogical);
p.addParameter('pathToPlotOutputDirRoot','~/Desktop/rgcDisplacementMapPlots',@ischar);

% parse
p.parse(varargin{:})


%% Setup
close all

% check or make a directory for output
if p.Results.savePlots
    if exist(fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName),'dir')==0
        mkdir(fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName));
    end
end

% Prepare the regular eccentricity support base
regularSupportPosDegRetina = ...
    0:p.Results.sampleResolutionDegrees:p.Results.maxModeledEccentricity;

% Define the cone and rgc density data to operate upon
coneDensityDataFileName = fullfile([getpref('rgcDisplacementMap','LocalDataPath') , '/Curcio_1990_JCompNeurol_HumanPhotoreceptorTopography/curcioRawConeDensity_',p.Results.subjectName,'.mat']);
rgcDensityDataFileName = fullfile([getpref('rgcDisplacementMap','LocalDataPath') , '/Curcio_1990_JCompNeurol_GanglionCellTopography/curcioRawRGCDensity_',p.Results.subjectName,'.mat']);

% Get the displacement map
[ ~, fitParams, meridianAngleSupport, rgcDisplacementEachMeridian, mRGC_cumulativeEachMeridian, mRF_cumulativeEachMeridian ] = ...
    makeDisplacementMap(...
    'sampleResolutionDegrees', p.Results.sampleResolutionDegrees, ...
    'maxModeledEccentricity', p.Results.maxModeledEccentricity, ...
    'meridianAngleResolutionDeg', p.Results.meridianAngleResolutionDeg, ...
    'displacementMapPixelsPerDeg', p.Results.displacementMapPixelsPerDeg, ...
    'coneDensityDataFileName', coneDensityDataFileName, ...
    'rgcDensityDataFileName', rgcDensityDataFileName, ...
    'verbose', p.Results.verbose);


%% Loop over the meridians
% Create cone density and RGC density polar maps
for mm = 1:length(meridianAngleSupport)
    
    % obtain cone density
    fitConeDensitySqDegRetina = getSplineFitToConeDensitySqDegRetina(meridianAngleSupport(mm));
    coneDensitySqDeg = zeroOpticDiscPoints(fitConeDensitySqDegRetina(regularSupportPosDegRetina),regularSupportPosDegRetina, meridianAngleSupport(mm));
    coneDensityEachMeridian(mm,:) = coneDensitySqDeg;
    
    % obtain the mRF density
    [ mRFDensitySqDeg, mRFtoConeDensityRatio ] = transformConeToMidgetRFDensity( coneDensitySqDeg, 'linkingFuncParams', fitParams(mm,1:2) );
    mRFDensityEachMeridian(mm,:) = mRFDensitySqDeg;
    mRFtoConeDensityEachMeridian(mm,:) = mRFtoConeDensityRatio;
    
    % obtain the RGC density
    fitRGCDensitySqDegRetina = getSplineFitToRGCDensitySqDegRetina(meridianAngleSupport(mm));
    rgcDensitySqDeg = zeroOpticDiscPoints(fitRGCDensitySqDegRetina(regularSupportPosDegRetina),regularSupportPosDegRetina, meridianAngleSupport(mm));
    rgcDensityEachMeridian(mm,:) = rgcDensitySqDeg;
    
    % obtain the mRGC density
    [ mRGCDensitySqDeg, midgetFraction ] = transformRGCToMidgetRGCDensityDacey( regularSupportPosDegRetina, rgcDensitySqDeg, 'linkingFuncParams', fitParams(mm,3:end) );
    mRGCDensityEachMeridian(mm,:) = mRGCDensitySqDeg;
    midgetFractionEachMeridian(mm,:) = midgetFraction;
    
end

% Define the image sample base for transform from polar coords
imRdim = (p.Results.maxModeledEccentricity * p.Results.displacementMapPixelsPerDeg * 2) -1;

% Show and save the maps
polarMapNameList = {...
    'coneDensityEachMeridian',...
    'rgcDensityEachMeridian',...
    'mRFDensityEachMeridian',...
    'mRFtoConeDensityEachMeridian',...
    'mRGCDensityEachMeridian',...
    'midgetFractionEachMeridian',...
    'rgcDisplacementEachMeridian',...
    'mRGC_cumulativeEachMeridian',...
    'mRF_cumulativeEachMeridian'
    };

% loop over the maps
for vv = 1:length(polarMapNameList)
    mapImage = feval('convertPolarMapToImageMap', eval(polarMapNameList{vv}), imRdim);
    figHandle = figure();
    climVals = [0,ceil(max(max(mapImage)))];
    tmp = strsplit(polarMapNameList{vv},'EachMeridian');
    titleString=tmp{1};
    displayRetinalImage(mapImage, climVals, p.Results.displacementMapPixelsPerDeg, p.Results.maxModeledEccentricity, titleString);
    if p.Results.savePlots
        fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,[titleString '.pdf']);
        saveas(figHandle,fileOutPath)
        close(figHandle);
    end
end

% warp some of the maps to cone space
warpMapNameList = {...
    'rgcDensityEachMeridian',...
    'mRGCDensityEachMeridian',...
    'mRGC_cumulativeEachMeridian',...
    };

% create a sample space for the warped map
eccenExtent = p.Results.maxModeledEccentricity - (1/p.Results.displacementMapPixelsPerDeg)/2;
smps = -eccenExtent+(1/p.Results.displacementMapPixelsPerDeg)/2:1/p.Results.displacementMapPixelsPerDeg:eccenExtent-(1/p.Results.displacementMapPixelsPerDeg)/2;
[sampleBaseX,sampleBaseY] = meshgrid(smps,smps);

% obtain the displacement map
displacementMapDeg = convertPolarMapToImageMap( rgcDisplacementEachMeridian, imRdim);

% loop over the maps
for vv = 1:length(warpMapNameList)
    mapImage = feval('convertPolarMapToImageMap', eval(warpMapNameList{vv}), imRdim);
    warpImage = applyDisplacementMap( mapImage, displacementMapDeg, sampleBaseX, sampleBaseY );
    smoothImage = fillAndSmoothMap(warpImage,sampleBaseX,sampleBaseY);
    figHandle = figure();
    climVals = [0,ceil(max(max(smoothImage)))];
    tmp = strsplit(warpMapNameList{vv},'EachMeridian');
    displayRetinalImage(smoothImage, climVals, p.Results.displacementMapPixelsPerDeg, p.Results.maxModeledEccentricity, ['warped ' tmp{1} ]);
    if p.Results.savePlots
        fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,['warped' titleString '.pdf']);
        saveas(figHandle,fileOutPath)
        close(figHandle);
    end
end % loop over maps to warp

% Make an image which is the difference between the warped mRGC_cumulative
% map and the mRF_cumulative
mapImageA = convertPolarMapToImageMap(mRF_cumulativeEachMeridian, imRdim);
mapImageB = applyDisplacementMap( ...
    convertPolarMapToImageMap(mRGC_cumulativeEachMeridian, imRdim), ...
    displacementMapDeg, sampleBaseX, sampleBaseY);
mapImageBminusA = mapImageB-mapImageA;
smoothImage = fillAndSmoothMap(mapImageBminusA,sampleBaseX,sampleBaseY);
figHandle = figure();
climVals = [0, 1e4];
titleString='mRGC_cumulative_warped_minus_mRF_cumulative';
displayRetinalImage(smoothImage, climVals, p.Results.displacementMapPixelsPerDeg, p.Results.maxModeledEccentricity, titleString);
if p.Results.savePlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,[titleString '.pdf']);
    saveas(figHandle,fileOutPath)
    close(figHandle);
end

mapImageAminusB = mapImageA-mapImageB;
smoothImage = fillAndSmoothMap(mapImageAminusB,sampleBaseX,sampleBaseY);
figHandle = figure();
climVals = [0, 1e4];
titleString='mRF_cumulative_minus_mRGC_cumulative_warped';
displayRetinalImage(smoothImage, climVals, p.Results.displacementMapPixelsPerDeg, p.Results.maxModeledEccentricity, titleString);
if p.Results.savePlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,[titleString '.pdf']);
    saveas(figHandle,fileOutPath)
    close(figHandle);
end

end % function


%% LOCAL FUNCTIONS

function vectorOut = zeroOpticDiscPoints(vectorIn, regularSupportPosDegRetina, polarAngle)
	opticDiscIndices = findOpticDiscPositions(regularSupportPosDegRetina, polarAngle);
    vectorOut = vectorIn;
    vectorOut(opticDiscIndices) = 0;
end
