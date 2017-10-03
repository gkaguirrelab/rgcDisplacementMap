function validation_ShowAllMaps(varargin)

% Transform cone density map to RGC density map


%% Parse input and define variables
p = inputParser; p.KeepUnmatched = true;

% Optional anaysis params
p.addParameter('sampleResolutionDegrees',0.01,@isnumeric);
p.addParameter('maxModeledEccentricity',30,@isnumeric);
p.addParameter('meridianAngleResolutionDeg',1,@isnumeric);
p.addParameter('displacementMapPixelsPerDeg',10,@isnumeric);
p.addParameter('pathToPlotOutputDir','~/Desktop/rgcDisplacementMapPlots',@ischar);

% Optional display and ouput params
p.addParameter('verbose',false,@islogical);
p.addParameter('savePlots',true,@islogical);

% parse
p.parse(varargin{:})


%% Setup
close all

% check or make a directory for output
if p.Results.savePlots
    if exist(p.Results.pathToPlotOutputDir,'dir')==0
        mkdir(p.Results.pathToPlotOutputDir);
    end
end

% Prepare the regular eccentricity support base
regularSupportPosDeg = ...
    0:p.Results.sampleResolutionDegrees:p.Results.maxModeledEccentricity;

% Get the displacement map
[ ~, fitParams, meridianAngles, rgcDisplacementEachMeridian, mRGC_cumulativeEachMeridian, mRF_cumulativeEachMeridian ] = ...
    makeDisplacementMap(...
    'sampleResolutionDegrees', p.Results.sampleResolutionDegrees, ...
    'maxModeledEccentricity', p.Results.maxModeledEccentricity, ...
    'meridianAngleResolutionDeg', p.Results.meridianAngleResolutionDeg, ...
    'displacementMapPixelsPerDeg', p.Results.displacementMapPixelsPerDeg);


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

% Define the image sample base for transform from polar coords
imRdim = p.Results.maxModeledEccentricity * p.Results.displacementMapPixelsPerDeg * 2;

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
    imagesc(mapImage, climVals);
    axis square
    set(gca,'TickLength',[0 0])
    tmp = strsplit(polarMapNameList{vv},'EachMeridian');
    titleString=tmp{1};
    c = colorbar;
    c.Label.String=[titleString];
    xlabel('Position [deg] nasal --> temporal');
    ylabel('Position [deg] inferior --> superior');
    numTicks=length(xticks);
    k=linspace(-1*p.Results.maxModeledEccentricity,p.Results.maxModeledEccentricity,numTicks+1);
    xticklabels(string(k(2:end)))
    yticklabels(string(k(2:end)))
    if p.Results.savePlots
        fileOutPath = fullfile(p.Results.pathToPlotOutputDir,[titleString '.pdf']);
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
smps = -eccenExtent:1/p.Results.displacementMapPixelsPerDeg:eccenExtent;
[sampleBaseX,sampleBaseY] = meshgrid(smps,smps);

% obtain the displacement map
displacementMapDeg = convertPolarMapToImageMap( rgcDisplacementEachMeridian, imRdim);

% loop over the maps
for vv = 1:length(warpMapNameList)
    mapImage = feval('convertPolarMapToImageMap', eval(warpMapNameList{vv}), imRdim);
    warpImage = applyDisplacementMap( mapImage, displacementMapDeg, sampleBaseX, sampleBaseY );
    figHandle = figure();
    climVals = [0,ceil(max(max(warpImage)))];
    imagesc(warpImage, climVals);
    axis square
    set(gca,'TickLength',[0 0])
    tmp = strsplit(warpMapNameList{vv},'EachMeridian');
    titleString=tmp{1};
    c = colorbar;
    c.Label.String=['warped ' titleString ];
    xlabel('Position [deg] nasal --> temporal');
    ylabel('Position [deg] inferior --> superior');
    numTicks=length(xticks);
    k=linspace(-1*p.Results.maxModeledEccentricity,p.Results.maxModeledEccentricity,numTicks+1);
    xticklabels(string(k(2:end)))
    yticklabels(string(k(2:end)))
    if p.Results.savePlots
        fileOutPath = fullfile(p.Results.pathToPlotOutputDir,['warped' titleString '.pdf']);
        saveas(figHandle,fileOutPath)
        close(figHandle);
    end
end % loop over maps to warp

% Make an image which is the difference between the mRF_cumulative and the
% warped mRGC_cumulative maps
mapImageA = convertPolarMapToImageMap(mRF_cumulativeEachMeridian, imRdim);
mapImageB = applyDisplacementMap( ...
    convertPolarMapToImageMap(mRGC_cumulativeEachMeridian, imRdim), ...
    displacementMapDeg, sampleBaseX, sampleBaseY);
mapImage = mapImageA - mapImageB;
figHandle = figure();
climVals = [floor(min(min(mapImage))),ceil(max(max(mapImage)))];
imagesc(mapImage, climVals);
axis square
set(gca,'TickLength',[0 0])
titleString='mRF_minus_mRGCwarped';
c = colorbar;
c.Label.String= titleString ;
xlabel('Position [deg] nasal --> temporal');
ylabel('Position [deg] inferior --> superior');
numTicks=length(xticks);
k=linspace(-1*p.Results.maxModeledEccentricity,p.Results.maxModeledEccentricity,numTicks+1);
xticklabels(string(k(2:end)))
yticklabels(string(k(2:end)))
if p.Results.savePlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDir,[titleString '.pdf']);
    saveas(figHandle,fileOutPath)
    close(figHandle);
end

end % function


function imageMap = convertPolarMapToImageMap(polarMap, imRdim)
maxMapDensity = max(polarMap(:));
imP=polarMap'./maxMapDensity;
imR = PolarToIm (imP, 0, 1, imRdim, imRdim);
imageMap = imrotate(imR .* maxMapDensity,-90);
end

