function [rgcDensitySqDeg, supportPosDeg] = loadRawRGCDensityByEccen(polarAngle, varargin)
% loadRawRGCDensityByEccen(polarAngle)
%
% This routine loads RGC density data from a file in standardized format.
% The file must be a .mat file that contains a single variable that it
% itself a structure, containing the fields:
%
%   support
%   inferior
%   superior
%   nasal
%   temporal
%   meta
%
% All fields except for meta contain a [1 x n] row-array of values. All
% these array fields must be of the same size.
%
% The meta field must contain entries for densityUnits and supportUnits,
% and these entries must be one of the approved types.
%
%
% Inputs:
%   polarAngle - The desired angle of the density function on the retinal field.
%                (0=nasal;90=superior;180=temporal;270=inferior)
%
% Outputs:
%   RGCDensitySqDeg - the density (counts per square degree) of RGCs at
%       each of the positions
%   supportPosDeg - the positions (in degrees of visual angle) from the
%       fovea at which the RGC density is defined
%
% Options:
%  densityDataFileName - The full path to the data file. The default value
%       assigned here is the average RGC density reported in Curcio 1990.


%% Parse input and define variables
p = inputParser;

% Required
p.addRequired('polarAngle',@isnumeric);

% Optional anaysis params
p.addParameter('densityDataFileName', ...
    fullfile([getpref('rgcDisplacementMap','LocalDataPath') , '/Curcio_1990_JCompNeurol_GanglionCellTopography/curcioRawRGCDensity_average.mat']), ...
    @ischar);
p.addParameter('cardinalMeridianAngles',[0 90 180 270],@isnumeric);
p.addParameter('cardinalMeridianNames',{'nasal' 'superior' 'temporal' 'inferior'},@iscell);

% parse
p.parse(polarAngle, varargin{:})

% Check the input
if sum(p.Results.cardinalMeridianAngles==polarAngle) ~= 1
    error('Please identify the polar angle for a cardinal meridian');
end

%% Load and check the density data
tmpLoader = load(p.Results.densityDataFileName);
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

%% Select the requested meridian and perform unit conversion
switch rawRGCDensity.meta.supportUnits
    case {'mm','MM','Mm'}
        % Convert mm to deg
        supportPosDeg = ...
            convert_mm_to_deg(rawRGCDensity.support);
    case {'deg','degrees'}
        % no conversion needed
        supportPosDeg = rawRGCDensity.support;
    otherwise
        error('The supportUnits of this raw file are not recognized');
end

% Place the requested meridian in a temporary variable
requestedMeridianName = p.Results.cardinalMeridianNames{find(p.Results.cardinalMeridianAngles == polarAngle)};
rawRGCDensityForSelectedMeridian = rawRGCDensity.(requestedMeridianName);

% Perform denisty unit conversion
switch rawRGCDensity.meta.densityUnits
    case 'counts/mm2'
        rgcDensitySqDeg = ...
            convert_mmSq_to_degSq(supportPosDeg, rawRGCDensityForSelectedMeridian );
    case 'counts/deg2'
        % no conversion needed
        rgcDensitySqDeg = rawRGCDensityForSelectedMeridian;
    otherwise
        error('The densityUnits of this raw file are not recognized');
end


end % function

