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
p.addParameter('sampleResolutionDegreesRetina',0.01,@isnumeric);
p.addParameter('maxModeledEccentricityDegreesRetina',30,@isnumeric);
p.addParameter('meridianAngleResolutionDeg',30,@isnumeric);
p.addParameter('displacementMapPixelsPerDegRetina',10,@isnumeric);
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
[ rgcDisplacementByMeridian, meridianAngleSupport, regularSupportPosDegRetina, opticDiscLocationByMeridian, mRF_RingCumulativeByMeridian, mRGC_RingCumulativeByMeridian, fitParamsByMeridian, ~, convergenceEccenDegRetinaByMeridian ] = ...
    createDisplacementModel(...
    'sampleResolutionDegreesRetina', p.Results.sampleResolutionDegreesRetina, ...
    'maxModeledEccentricityDegreesRetina', p.Results.maxModeledEccentricityDegreesRetina, ...
    'meridianAngleResolutionDeg', p.Results.meridianAngleResolutionDeg, ...
    'displacementMapPixelsPerDegRetina', p.Results.displacementMapPixelsPerDegRetina, ...
    'coneDensityDataFileName', coneDensityDataFileName, ...
    'rgcDensityDataFileName', rgcDensityDataFileName, ...
    'verbose', p.Results.verbose);


%% Loop over the meridians
% Create cone density and RGC density polar maps
for mm = 1:length(meridianAngleSupport)
    
    % obtain cone density
    fitConeDensitySqDegRetina = getSplineFitToConeDensitySqDegRetina(meridianAngleSupport(mm));
    coneDensitySqDeg = zeroOpticDiscPoints(fitConeDensitySqDegRetina(regularSupportPosDegRetina),regularSupportPosDegRetina, meridianAngleSupport(mm));
    coneDensityByMeridian(mm,:) = coneDensitySqDeg;
    
    % obtain the mRF density
    [ mRFDensitySqDeg, mRFtoConeDensityRatio ] = transformConeToMidgetRFDensity( coneDensitySqDeg, 'linkingFuncParams', fitParamsByMeridian(mm,1:2) );
    mRFDensityByMeridian(mm,:) = mRFDensitySqDeg;
    mRFtoConeDensityByMeridian(mm,:) = mRFtoConeDensityRatio;
    
    % obtain the RGC density
    fitRGCDensitySqDegRetina = getSplineFitToRGCDensitySqDegRetina(meridianAngleSupport(mm));
    rgcDensitySqDeg = zeroOpticDiscPoints(fitRGCDensitySqDegRetina(regularSupportPosDegRetina),regularSupportPosDegRetina, meridianAngleSupport(mm));
    rgcDensityByMeridian(mm,:) = rgcDensitySqDeg;
    
    % obtain the mRGC density
    [ mRGCDensitySqDeg, midgetFraction ] = transformRGCToMidgetRGCDensityDacey( regularSupportPosDegRetina, rgcDensitySqDeg, 'linkingFuncParams', fitParamsByMeridian(mm,3:end) );
    mRGCDensityByMeridian(mm,:) = mRGCDensitySqDeg;
    midgetFractionByMeridian(mm,:) = midgetFraction;
    
end

% Define the image sample base for transform from polar coords
imRdim = (p.Results.maxModeledEccentricityDegreesRetina * p.Results.displacementMapPixelsPerDegRetina * 2) -1;

% Show and save the maps
polarMapNameList = {'rgcDisplacementByMeridian',...
    };

% Obtain the boundary of the optic disc in image space. We have to reverse
% the boundary array for the y-axis to match our image convention
opticDiscLocationsImage = convertPolarMapToImageMap(opticDiscLocationByMeridian, 'imRdim', imRdim);
opticDiscBoundary = bwboundaries(imbinarize(opticDiscLocationsImage),'noholes');
opticDiscBoundaryArray = opticDiscBoundary{1};
opticDiscBoundaryArray(:,1)=imRdim(1)-opticDiscBoundaryArray(:,1);
dashIndices=1:10:size(opticDiscBoundaryArray,1);

% Obtain the boundary of the end of the displacement zone. 
% make matrix that is NaNs size of length(meridianAngleSupport) by
% length(regularSupportPosDegRetina)
dispZoneLocationByMeridian = zeros(length(meridianAngleSupport), length(regularSupportPosDegRetina));
for ii  = 1:length(convergenceEccenDegRetinaByMeridian)
    dispZoneEndIndex = find(regularSupportPosDegRetina == convergenceEccenDegRetinaByMeridian(ii));
    dispZoneLocationByMeridian(ii,1:dispZoneEndIndex) = 1;
end
dispZoneLocationsImage = convertPolarMapToImageMap(dispZoneLocationByMeridian, 'imRdim', imRdim);
dispZoneBoundary = bwboundaries(imbinarize(dispZoneLocationsImage),'noholes');
dispZoneBoundaryArray = dispZoneBoundary{1};
dispZoneBoundaryArray(:,1)=imRdim(1)-dispZoneBoundaryArray(:,1);
dispZoneDashIndices=1:10:size(dispZoneBoundaryArray,1);


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
    displayRetinalImage(mapImage, climVals, p.Results.displacementMapPixelsPerDegRetina, p.Results.maxModeledEccentricityDegreesRetina, titleString);
    hold on
    if median(arrayfun(@(x) mapImage(opticDiscBoundaryArray(x,1),opticDiscBoundaryArray(x,2)),1:1:size(opticDiscBoundaryArray,1))) > max(climVals)/2
        opticDiscBorderColor = [0.25 0.25 0.25];
    else
        opticDiscBorderColor = [0.75 0.75 0.75];
        
    end
    dispZoneBorderColor = 'r';
    plot(opticDiscBoundaryArray(dashIndices,2),opticDiscBoundaryArray(dashIndices,1),'.','Color',opticDiscBorderColor);
    plot(dispZoneBoundaryArray(dispZoneDashIndices,2),dispZoneBoundaryArray(dispZoneDashIndices,1),'.','Color',dispZoneBorderColor);
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



end % function


%% LOCAL FUNCTIONS

function vectorOut = zeroOpticDiscPoints(vectorIn, regularSupportPosDegRetina, polarAngle)
opticDiscIndices = findOpticDiscPositions(regularSupportPosDegRetina, polarAngle);
vectorOut = vectorIn;
vectorOut(opticDiscIndices) = 0;
end
