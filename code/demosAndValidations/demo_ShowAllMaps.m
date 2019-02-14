function demo_ShowAllMaps(varargin)
% Displays and saves maps of retinal cell density and RGC displacement
%
% Description:
%	This validation function runs the displacement model for each of many
%	meridian values and then converts the model output to images which are
%	displayed and optionally saved to disk. The analysis parameters are
%   explicitly defined here and then passed to the main,
%   createDisplacementModel routine.
%
%   The "subjectName" corresponds to one of the Curcio 1990 datasets
%   present within the data directory of this toolbox. Typical options
%   include:
%       reportedAverage   - The values reported in the Curcio 1990 papers
%       computedAverage   - Our derivation of average values from the 
%                           Curcio data
%       29986A            - Data from subject 29986, averaged over both
%                           eyes
%


%% Parse input and define variables
p = inputParser; p.KeepUnmatched = true;

% Optional anaysis params
p.addParameter('sampleResolutionDegVisual',0.01,@isnumeric);
p.addParameter('maxModeledEccentricityDegVisual',30,@isnumeric);
p.addParameter('meridianAngleResolutionDeg',1,@isnumeric);
p.addParameter('displacementMapPixelsPerDegVisual',10,@isnumeric);
p.addParameter('subjectName', 'computedAverage', @ischar);

% Optional display and ouput params
p.addParameter('verbose',true,@islogical);
p.addParameter('saveRasterPlots',true,@islogical);
p.addParameter('saveVectorPlots',false,@islogical);
p.addParameter('pathToPlotOutputDirRoot','~/Desktop/rgcDisplacementMapPlots',@ischar);

% parse
p.parse(varargin{:})


%% Setup
close all

% check or make a directory for output
if p.Results.saveRasterPlots || p.Results.saveVectorPlots
    if exist(fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName),'dir')==0
        mkdir(fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName));
    end
end


% Define the cone and rgc density data to operate upon
coneDensityDataFileName = fullfile([getpref('rgcDisplacementMap','LocalDataPath') , '/Curcio_1990_JCompNeurol_HumanPhotoreceptorTopography/curcioRawConeDensity_',p.Results.subjectName,'.mat']);
rgcDensityDataFileName = fullfile([getpref('rgcDisplacementMap','LocalDataPath') , '/Curcio_1990_JCompNeurol_GanglionCellTopography/curcioRawRGCDensity_',p.Results.subjectName,'.mat']);

% Create the displacement model
[ rgcDisplacementByMeridian, meridianAngleSupport, regularSupportPosDegVisual, opticDiscLocationByMeridian, mRF_RingCumulativeByMeridian, mRGC_RingCumulativeByMeridian, fitParamsByMeridian, ~, ~ ] = ...
    createDisplacementModel(...
    'sampleResolutionDegVisual', p.Results.sampleResolutionDegVisual, ...
    'maxModeledEccentricityDegVisual', p.Results.maxModeledEccentricityDegVisual, ...
    'meridianAngleResolutionDeg', p.Results.meridianAngleResolutionDeg, ...
    'displacementMapPixelsPerDegVisual', p.Results.displacementMapPixelsPerDegVisual, ...
    'coneDensityDataFileName', coneDensityDataFileName, ...
    'rgcDensityDataFileName', rgcDensityDataFileName, ...
    'verbose', p.Results.verbose);


%% Loop over the meridians
% Create cone density and RGC density polar maps
for mm = 1:length(meridianAngleSupport)
    
    % obtain cone density
    fitConeDensitySqDegVisual = getSplineFitToConeDensitySqDegVisual(meridianAngleSupport(mm));
    coneDensitySqDeg = zeroOpticDiscPoints(fitConeDensitySqDegVisual(regularSupportPosDegVisual),regularSupportPosDegVisual, meridianAngleSupport(mm));
    coneDensityByMeridian(mm,:) = coneDensitySqDeg;
    
    % obtain the mRF density
    [ mRFDensitySqDeg, mRFtoConeDensityRatio ] = transformConeToMidgetRFDensity( coneDensitySqDeg, 'linkingFuncParams', fitParamsByMeridian(mm,1:2) );
    mRFDensityByMeridian(mm,:) = mRFDensitySqDeg;
    mRFtoConeDensityByMeridian(mm,:) = mRFtoConeDensityRatio;
    
    % obtain the RGC density
    fitRGCDensitySqDegVisual = getSplineFitToRGCDensitySqDegVisual(meridianAngleSupport(mm));
    rgcDensitySqDeg = zeroOpticDiscPoints(fitRGCDensitySqDegVisual(regularSupportPosDegVisual),regularSupportPosDegVisual, meridianAngleSupport(mm));
    rgcDensityByMeridian(mm,:) = rgcDensitySqDeg;
    
    % obtain the mRGC density
    [ mRGCDensitySqDeg, midgetFraction ] = transformRGCToMidgetRGCDensityDacey( regularSupportPosDegVisual, rgcDensitySqDeg, 'linkingFuncParams', fitParamsByMeridian(mm,3:end) );
    mRGCDensityByMeridian(mm,:) = mRGCDensitySqDeg;
    midgetFractionByMeridian(mm,:) = midgetFraction;
    
end

% Define the image sample base for transform from polar coords
imRdim = (p.Results.maxModeledEccentricityDegVisual * p.Results.displacementMapPixelsPerDegVisual * 2) -1;

% Show and save the maps
polarMapNameList = {...
    'coneDensityByMeridian',...
    'rgcDensityByMeridian',...
    'mRFDensityByMeridian',...
    'mRFtoConeDensityByMeridian',...
    'mRGCDensityByMeridian',...
    'midgetFractionByMeridian',...
    'rgcDisplacementByMeridian',...
    'mRGC_RingCumulativeByMeridian',...
    'mRF_RingCumulativeByMeridian' ...
    };

% Obtain the boundary of the optic disc in image space. We have to reverse
% the boundary array for the y-axis to match our image convention
opticDiscLocationsImage = convertPolarMapToImageMap(opticDiscLocationByMeridian, 'imRdim', imRdim);
opticDiscBoundary = bwboundaries(imbinarize(opticDiscLocationsImage),'noholes');
opticDiscBoundaryArray = opticDiscBoundary{1};
opticDiscBoundaryArray(:,1)=imRdim(1)-opticDiscBoundaryArray(:,1);
dashIndices=1:10:size(opticDiscBoundaryArray,1);

% loop over the maps
for vv = 1:length(polarMapNameList)
    mapImage = feval('convertPolarMapToImageMap', eval(polarMapNameList{vv}), 'imRdim', imRdim);
    figHandle = figure();
    figHandle.Renderer='Painters';
    climVals = [0,ceil(max(max(mapImage)))];
    tmp = strsplit(polarMapNameList{vv},'ByMeridian');
    titleString=tmp{1};
    % Handle the special case of the clim for mRFtoConeDensityByMeridian
    if strcmp(titleString,'mRFtoConeDensity')
        climVals(2)=2;
    end
    displayRetinalImage(mapImage, climVals, p.Results.displacementMapPixelsPerDegVisual, p.Results.maxModeledEccentricityDegVisual, titleString);
    hold on
    if median(arrayfun(@(x) mapImage(opticDiscBoundaryArray(x,1),opticDiscBoundaryArray(x,2)),1:1:size(opticDiscBoundaryArray,1))) > max(climVals)/2
        opticDiscBorderColor = [0.25 0.25 0.25];
    else
        opticDiscBorderColor = [0.75 0.75 0.75];
        
    end
    plot(opticDiscBoundaryArray(dashIndices,2),opticDiscBoundaryArray(dashIndices,1),'.','Color',opticDiscBorderColor);
    if p.Results.saveRasterPlots
        % save a rasterized figure without the optic disc marker
        fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,[titleString '.png']);
        print(fileOutPath,'-dtiffn','-r300');
    end
    if p.Results.saveVectorPlots
        % save a vector figure with the optic disc marker
        fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,[titleString '.pdf']);
        saveas(figHandle,fileOutPath)
    end
    close(figHandle);
end

% warp some of the maps to cone space
warpMapNameList = {...
    'rgcDensityByMeridian',...
    'mRGCDensityByMeridian',...
    'mRGC_RingCumulativeByMeridian',...
    };

% create a sample space for the warped map
eccenExtent = p.Results.maxModeledEccentricityDegVisual - (1/p.Results.displacementMapPixelsPerDegVisual)/2;
smps = -eccenExtent+(1/p.Results.displacementMapPixelsPerDegVisual)/2:1/p.Results.displacementMapPixelsPerDegVisual:eccenExtent-(1/p.Results.displacementMapPixelsPerDegVisual)/2;
[sampleBaseX,sampleBaseY] = meshgrid(smps,smps);

% obtain the displacement map
displacementMapDeg = convertPolarMapToImageMap( rgcDisplacementByMeridian, 'imRdim', imRdim);

% loop over the maps
for vv = 1:length(warpMapNameList)
    mapImage = feval('convertPolarMapToImageMap', eval(warpMapNameList{vv}), 'imRdim', imRdim);
    warpImage = applyDisplacementMap( mapImage, displacementMapDeg, sampleBaseX, sampleBaseY );
    smoothImage = fillAndSmoothMap(warpImage,sampleBaseX,sampleBaseY);
    figHandle = figure();
    figHandle.Renderer='Painters';
    climVals = [0,ceil(max(max(smoothImage)))];
    tmp = strsplit(warpMapNameList{vv},'ByMeridian');
    titleString=tmp{1};
    displayRetinalImage(smoothImage, climVals, p.Results.displacementMapPixelsPerDegVisual, p.Results.maxModeledEccentricityDegVisual, ['warped ' tmp{1} ]);
    hold on
    if median(arrayfun(@(x) smoothImage(opticDiscBoundaryArray(x,1),opticDiscBoundaryArray(x,2)),1:1:size(opticDiscBoundaryArray,1))) > max(climVals)/2
        opticDiscBorderColor = [0.25 0.25 0.25];
    else
        opticDiscBorderColor = [0.75 0.75 0.75];
        
    end
    plot(opticDiscBoundaryArray(dashIndices,2),opticDiscBoundaryArray(dashIndices,1),'.','Color',opticDiscBorderColor);
    if p.Results.saveRasterPlots
        % save a rasterized figure without the optic disc marker
        fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,[titleString '_warped.png']);
        print(fileOutPath,'-dtiffn','-r300');
    end
    if p.Results.saveVectorPlots
        % save a vector figure with the optic disc marker
        fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,[titleString '_warped.pdf']);
        saveas(figHandle,fileOutPath)
    end
    close(figHandle);
end % loop over maps to warp

% Make an image which is the difference between the warped mRGC_cumulative
% map and the mRF_cumulative
mapImageA = convertPolarMapToImageMap(mRF_RingCumulativeByMeridian, 'imRdim', imRdim);
mapImageB = applyDisplacementMap( ...
    convertPolarMapToImageMap(mRGC_RingCumulativeByMeridian, 'imRdim', imRdim), ...
    displacementMapDeg, sampleBaseX, sampleBaseY);
mapImageBminusA = mapImageB-mapImageA;
smoothImage = fillAndSmoothMap(mapImageBminusA,sampleBaseX,sampleBaseY);
figHandle = figure();
figHandle.Renderer='Painters';
climVals = [0, 1e4];
titleString='mRGC_cumulative_warped_minus_mRF_cumulative';
displayRetinalImage(smoothImage, climVals, p.Results.displacementMapPixelsPerDegVisual, p.Results.maxModeledEccentricityDegVisual, titleString);
hold on
if median(arrayfun(@(x) smoothImage(opticDiscBoundaryArray(x,1),opticDiscBoundaryArray(x,2)),1:1:size(opticDiscBoundaryArray,1))) > max(climVals)/2
    opticDiscBorderColor = [0.25 0.25 0.25];
else
    opticDiscBorderColor = [0.75 0.75 0.75];
    
end
plot(opticDiscBoundaryArray(dashIndices,2),opticDiscBoundaryArray(dashIndices,1),'.','Color',opticDiscBorderColor);
if p.Results.saveRasterPlots
    % save a rasterized figure without the optic disc marker
    fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,[titleString '.png']);
    print(fileOutPath,'-dtiffn','-r300');
end
if p.Results.saveVectorPlots
    % save a vector figure with the optic disc marker
    fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,[titleString '.pdf']);
    saveas(figHandle,fileOutPath)
end
close(figHandle);

mapImageAminusB = mapImageA-mapImageB;
smoothImage = fillAndSmoothMap(mapImageAminusB,sampleBaseX,sampleBaseY);
figHandle = figure();
figHandle.Renderer='Painters';
climVals = [0, 1e4];
titleString='mRF_cumulative_minus_mRGC_cumulative_warped';
displayRetinalImage(smoothImage, climVals, p.Results.displacementMapPixelsPerDegVisual, p.Results.maxModeledEccentricityDegVisual, titleString);
hold on
if median(arrayfun(@(x) smoothImage(opticDiscBoundaryArray(x,1),opticDiscBoundaryArray(x,2)),1:1:size(opticDiscBoundaryArray,1))) > max(climVals)/2
    opticDiscBorderColor = [0.25 0.25 0.25];
else
    opticDiscBorderColor = [0.75 0.75 0.75];
    
end
plot(opticDiscBoundaryArray(dashIndices,2),opticDiscBoundaryArray(dashIndices,1),'.','Color',opticDiscBorderColor);
if p.Results.saveRasterPlots
    % save a rasterized figure without the optic disc marker
    fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,[titleString '.png']);
    print(fileOutPath,'-dtiffn','-r300');
end
if p.Results.saveVectorPlots
    % save a vector figure with the optic disc marker
    fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,[titleString '.pdf']);
    saveas(figHandle,fileOutPath)
end
close(figHandle);


end % function


%% LOCAL FUNCTIONS

function vectorOut = zeroOpticDiscPoints(vectorIn, regularSupportPosDegVisual, polarAngle)
opticDiscIndices = findOpticDiscPositions(regularSupportPosDegVisual, polarAngle);
vectorOut = vectorIn;
vectorOut(opticDiscIndices) = 0;
end
