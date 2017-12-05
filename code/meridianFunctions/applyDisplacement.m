function valuesAtCones = applyDisplacement(valuesAtRGCSoma, displacementDegRetina, supportPosDegRetina, varargin)
% Transforms a vector values at RGC locations to values at cone locations
%
% Description
%   Implements the displacement of values from RGC to cone locations.
%
% Inputs:
%   valuesAtRGCSoma       - A 1 x p vector of values at each of the p
%                           eccentricity locations.
%   displacementDegRetina - A 1 x p vector, where p is the number of
%                           eccentricity positions modeled. The values are
%                           the displacement in retinal degrees of the RGC
%                           soma away from the fovea.
%   supportPosDegRetina   - A 1 x p vector that contains the locations in
%                           eccentricity in retinal degrees for each of
%                           the valuesAtRGCSoma and displacementDegRetina
%
% Optional key/value pairs:
%  'sampleResolutionDegreesRetina' - The calculations are performed across
%                           a regular sampling of eccentricity. This param
%                           defines sample resolution The sample resolution
%                           must be sufficient fine so that the cumulative
%                           is an accurate estimate of the integral.
%
% Outputs:
%
%   valuesAtCones         - A 1 x p vector of values at each of the p
%                           eccentricity locations.
%

%% Parse input and define variables
p = inputParser;

% required input
p.addRequired('valuesAtRGCSoma',@isnumeric);
p.addRequired('displacementDegRetina',@isnumeric);
p.addRequired('supportPosDegRetina',@isnumeric);

% Optional anaysis params
p.addParameter('sampleResolutionDegreesRetina',0.01,@isnumeric);

% parse
p.parse(valuesAtRGCSoma, displacementDegRetina, supportPosDegRetina, varargin{:})

% Implement the transformation

% We first identfy the destination (in eccentricities) at each of the cone
% locations for each of the RGC locations (rounded to the closest support
% position)
destination=round(supportPosDegRetina-displacementDegRetina,ceil(-log10(p.Results.sampleResolutionDegreesRetina)));
uniqueDestinations = unique(destination);
destinationIdx=arrayfun(@(x) find((abs(supportPosDegRetina-x))<1e-6,1), uniqueDestinations, 'UniformOutput', false);
validDestinationIndices = cellfun(@(x) ~isempty(x), destinationIdx);
valuesAtConesAtValidDestinations = arrayfun(@(x) sum(valuesAtRGCSoma(destination==x)), uniqueDestinations(validDestinationIndices));
valuesAtCones=zeros(size(valuesAtRGCSoma));
valuesAtCones(unique(cell2mat(destinationIdx(validDestinationIndices))))=valuesAtConesAtValidDestinations;

end