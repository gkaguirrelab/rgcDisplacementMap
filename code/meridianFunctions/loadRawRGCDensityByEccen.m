function [rgcDensitySqDegVisual, supportPosDegVisual] = loadRawRGCDensityByEccen(polarAngle, varargin)
% Load RGC density data from standard file and return in visual degrees
%
% Syntax:
%  [rgcDensitySqDegVisual, supportPosDegField] = loadRawRGCDensityByEccen(polarAngle, varargin)
%
% Description:
%   This routine loads RGC density data from a file in standardized format.
%   The file must be a .mat file that contains a single variable that is
%   itself a structure, containing the fields:
%
%       support
%       inferior
%       superior
%       nasal
%       temporal
%       meta
%
%   All fields except for meta contain a [1 x n] row-array of values. All
%   these array fields must be of the same size.
%
%   The meta field must contain entries for densityUnits and supportUnits,
%   and these entries must be one of the approved types.
%
%   The output of the function is in units of count per square visual field
%   degree and in position of visual field
%
% Inputs:
%   polarAngle            - The desired angle of the density function on
%                           the retinal field. (0=nasal; 90=superior;
%                           180=temporal; 270=inferior)
%
% Optional key/value pairs:
%  rgcDensityDataFileName - The full path to the data file. The default
%                           value assigned here is the average RGC density
%                           reported in Curcio 1990.
%
% Outputs:
%   rgcDensitySqDegRetina - The density (counts per square degree of visual
%                           angle) of RGCs at each of the positions
%   supportPosDegRetina   - The positions (in degrees of visual angle)
%                           from the fovea at which the RGC density is
%                           defined
%
% Examples:
%{
    [rgcDensitySqDegVisual, supportPosDegVisual] = loadRawRGCDensityByEccen(180);
    plot(supportPosDegVisual, rgcDensitySqDegVisual,'*r')
    xlabel('Eccentricty [visual degrees]')
    ylabel('Total RGC density [counts / sq degree visual]');
    legend('temporal');
%}


%% Parse input and define variables
p = inputParser;

% Required
p.addRequired('polarAngle',@isnumeric);

% Optional anaysis params
p.addParameter('rgcDensityDataFileName', [],@(x)(isempty(x) | ischar(x)));
p.addParameter('cardinalMeridianAngles',[0 90 180 270],@isnumeric);
p.addParameter('cardinalMeridianNames',{'nasal' 'superior' 'temporal' 'inferior'},@iscell);

% parse
p.parse(polarAngle, varargin{:})

% Check the input
if sum(p.Results.cardinalMeridianAngles==polarAngle) ~= 1
    error('Please identify the polar angle for a cardinal meridian');
end

%% Set the default for the rgcDensityFileName if not defined
if isempty(p.Results.rgcDensityDataFileName)
    rgcDensityDataFileName = ...
        fullfile([getpref('rgcDisplacementMap','LocalDataPath') , '/Curcio_1990_JCompNeurol_GanglionCellTopography/curcioRawRGCDensity_computedAverage.mat']);
else
	rgcDensityDataFileName = p.Results.rgcDensityDataFileName;
end

%% Load and check the density data
tmpLoader = load(rgcDensityDataFileName);
fieldNames = fields(tmpLoader);
if length(fieldNames) ~= 1
    error('Please specify a raw datafile that contains a single structure variable');
end
rawRGCDensity = tmpLoader.(fieldNames{1});

% Check that the units meta fields are defined
if ~isfield(rawRGCDensity.meta, 'supportUnits')
    error('The raw data file lacks the field .meta.supportUnits');
end
if ~isfield(rawRGCDensity.meta, 'densityUnits')
    error('The raw data file lacks the field .meta.densityUnits');
end

% What is the meridian we are using?
requestedMeridianName = p.Results.cardinalMeridianNames{find(p.Results.cardinalMeridianAngles == polarAngle)};


%% Select the requested meridian and perform unit conversion
switch rawRGCDensity.meta.supportUnits
    case {'mm','MM','Mm'}
        if isstruct(rawRGCDensity.support)
            supportPosDegVisual = convert_mmRetina_to_degVisual(rawRGCDensity.support.(requestedMeridianName), polarAngle);
        else
            supportPosDegVisual = convert_mmRetina_to_degVisual(rawRGCDensity.support, polarAngle);
        end
    case {'visualDegree'}
        % no conversion needed
        if isstruct(rawRGCDensity.support)
            supportPosDegVisual = rawRGCDensity.support.(requestedMeridianName);
        else
            
            supportPosDegVisual = rawRGCDensity.support;
        end
    otherwise
        error('The supportUnits of this raw file are not recognized');
end

% Place the requested meridian in a temporary variable
rawRGCDensityForSelectedMeridian = rawRGCDensity.(requestedMeridianName);

% Perform density unit conversion
switch rawRGCDensity.meta.densityUnits
    case 'counts/mm2'
        if isstruct(rawRGCDensity.support)
            degSqVisualPerMmSqRetinaRelativeToVisualAxis = calc_degSqVisual_per_mmSqRetina(rawRGCDensity.support.(requestedMeridianName), polarAngle);
        else
            degSqVisualPerMmSqRetinaRelativeToVisualAxis = calc_degSqVisual_per_mmSqRetina(rawRGCDensity.support, polarAngle);
        end
        rgcDensitySqDegVisual = rawRGCDensityForSelectedMeridian ./ degSqVisualPerMmSqRetinaRelativeToVisualAxis;
    case 'counts/visualDeg2'
        rgcDensitySqDegVisual = rawRGCDensityForSelectedMeridian;
    otherwise
        error('The densityUnits of this raw file are not recognized');
end


end % function

