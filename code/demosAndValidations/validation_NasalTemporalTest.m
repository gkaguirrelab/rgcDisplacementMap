function validation_NasalTemporalTest(varargin)

% Steps through each routine and checks the returned values for the nasal
% and temporal meridians against published values


%% Parse input and define variables
p = inputParser;

% Optional anaysis params
p.addParameter('sampleResolutionDegrees',0.01,@isnumeric);
p.addParameter('maxModeledEccentricity',30,@isnumeric);
p.addParameter('nasalMeridianAngle',0,@isnumeric);
p.addParameter('temporalMeridianAngle',180,@isnumeric);
p.addParameter('meridianAngleResolutionDeg',15,@isnumeric);
p.addParameter('displacementMapPixelsPerDeg',10,@isnumeric);

% Optional display params
p.addParameter('verbose',true,@islogical);

% parse
p.parse(varargin{:})

close all

%% Setup
% Prepare the regular eccentricity support base
regularSupportPosDeg = ...
    0:p.Results.sampleResolutionDegrees:p.Results.maxModeledEccentricity;


%% Check cone density values
fprintf('Checking the function getCurcioConeDensityByEccen:\n');

% Test the nasal meridian
[rgcDensitySqDeg, supportPosDeg] = getCurcioConeDensityByEccen(p.Results.nasalMeridianAngle);
checkEccenMm = 1; checkDensityMmSq = 20962; % These values copied from Curcio data table
checkEccenDeg = convert_mm_to_deg(checkEccenMm);
checkDensityMmDeg = convert_mmSq_to_degSq(checkEccenDeg, checkDensityMmSq );
Idx = find(supportPosDeg==checkEccenDeg,1);
if checkDensityMmDeg == rgcDensitySqDeg(Idx)
    fprintf('\tNasal cone density value matches Curcio table\n');
else
    error('Nasal cone density value DOES NOT MATCH Curcio table\n');
end

% Test the temporal meridian
[rgcDensitySqDeg, supportPosDeg] = getCurcioConeDensityByEccen(p.Results.temporalMeridianAngle);
checkEccenMm = 1; checkDensityMmSq = 19699; % These values copied from Curcio data table
checkEccenDeg = convert_mm_to_deg(checkEccenMm);
checkDensityMmDeg = convert_mmSq_to_degSq(checkEccenDeg, checkDensityMmSq );
Idx = find(supportPosDeg==checkEccenDeg,1);
if checkDensityMmDeg == rgcDensitySqDeg(Idx)
    fprintf('\tTemporal cone density value matches Curcio table\n');
else
    error('Temporal cone density value DOES NOT MATCH Curcio table\n');
end


%% Check RGC density values
fprintf('Checking the function getCurcioRGCDensityByEccen:\n');

% Test the nasal meridian
[rgcDensitySqDeg, supportPosDeg] = getCurcioRGCDensityByEccen(p.Results.nasalMeridianAngle);
checkEccenMm = 1.0085; checkDensityMmSq = 31544.2; % These values copied from Curcio data table
checkEccenDeg = convert_mm_to_deg(checkEccenMm);
checkDensityMmDeg = convert_mmSq_to_degSq(checkEccenDeg, checkDensityMmSq );
Idx = find(supportPosDeg==checkEccenDeg,1);
if checkDensityMmDeg == rgcDensitySqDeg(Idx)
    fprintf('\tNasal RGC density value matches Curcio table\n');
else
    fprintf('Nasal RGC density value DOES NOT MATCH Curcio table\n');
end

% Test the temporal meridian
[rgcDensitySqDeg, supportPosDeg] = getCurcioRGCDensityByEccen(p.Results.temporalMeridianAngle);
checkEccenMm = 1.0085; checkDensityMmSq = 26921.8; % These values copied from Curcio data table
checkEccenDeg = convert_mm_to_deg(checkEccenMm);
checkDensityMmDeg = convert_mmSq_to_degSq(checkEccenDeg, checkDensityMmSq );
Idx = find(supportPosDeg==checkEccenDeg,1);
if checkDensityMmDeg == rgcDensitySqDeg(Idx)
    fprintf('\tTemporal RGC density value matches Curcio table\n');
else
    fprintf('Temporal RGC density value DOES NOT MATCH Curcio table\n');
end

%% Check Watson equation for mRF density
fprintf('Checking the function calcWatsonMidgetRFDensityByEccen:\n');
% It should be the case that the density at 50° for the true nasal meridian
% is greater than that for the temporal meridian
if calcWatsonMidgetRFDensityByEccen(50, p.Results.nasalMeridianAngle) > ...
        calcWatsonMidgetRFDensityByEccen(50, p.Results.temporalMeridianAngle)
    fprintf('\tTemporal and Nasal mRF density values are appropriate.\n');
else
    fprintf('\tTemporal and Nasal mRF density values are swapped.\n');
end


end % validation function