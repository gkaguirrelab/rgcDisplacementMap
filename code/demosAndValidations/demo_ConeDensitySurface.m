function demo_ConeDensitySurface( varargin )
% Create a curved, surface map of cone density
%
% Description:
%   We do this.
%





%% Parse input and define variables
p = inputParser; p.KeepUnmatched = true;

% Optional anaysis params
p.addParameter('sampleResolutionDegreesRetina',0.01,@isnumeric);
p.addParameter('maxModeledEccentricityDegreesRetina',30,@isnumeric);
p.addParameter('meridianAngleResolutionDeg',1,@isnumeric);
p.addParameter('displacementMapPixelsPerDegRetina',10,@isnumeric);
p.addParameter('subjectName', 'computedAverage', @ischar);

% Optional display and ouput params
p.addParameter('verbose',true,@islogical);
p.addParameter('saveRasterPlots',true,@islogical);
p.addParameter('saveVectorPlots',false,@islogical);
p.addParameter('pathToPlotOutputDirRoot','~/Desktop/rgcDisplacementMapPlots',@ischar);

% parse
p.parse(varargin{:})

close all

regularSupportPosDegRetina = 0:1:50;
meridianAngleSupport = 0:15:359;

%% Loop over the meridians
% Create cone density and RGC density polar maps
for mm = 1:length(meridianAngleSupport)
    fprintf('.');
    % obtain cone density
    fitConeDensitySqDegRetina = getSplineFitToConeDensitySqDegRetina(meridianAngleSupport(mm));
    coneDensitySqDeg = zeroOpticDiscPoints(fitConeDensitySqDegRetina(regularSupportPosDegRetina),regularSupportPosDegRetina, meridianAngleSupport(mm));
    coneDensityByMeridian(mm,:) = coneDensitySqDeg;
        
    % obtain the RGC density
    fitRGCDensitySqDegRetina = getSplineFitToRGCDensitySqDegRetina(meridianAngleSupport(mm));
    rgcDensitySqDeg = zeroOpticDiscPoints(fitRGCDensitySqDegRetina(regularSupportPosDegRetina),regularSupportPosDegRetina, meridianAngleSupport(mm));
    rgcDensityByMeridian(mm,:) = rgcDensitySqDeg;
        
end
fprintf('\n');

[Xout, Yout, Zout, Cmap] = sphere3d(fliplr(flipud(rgcDensityByMeridian))',0,deg2rad(23*15),deg2rad(-90),deg2rad(-60),10,.25,'surf',1e-6);
set(findobj(gcf, 'type', 'surface'),'edgecolor','none');

end




%% LOCAL FUNCTIONS

function vectorOut = zeroOpticDiscPoints(vectorIn, regularSupportPosDegRetina, polarAngle)
opticDiscIndices = findOpticDiscPositions(regularSupportPosDegRetina, polarAngle);
vectorOut = vectorIn;
vectorOut(opticDiscIndices) = 0;
end
