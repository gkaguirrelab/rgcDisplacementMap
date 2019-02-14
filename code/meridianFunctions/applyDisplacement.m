function valuesAtCones = applyDisplacement(valuesAtRGCSoma, displacementDegVisual, supportPosDegVisual, varargin)
% Transforms a vector values at RGC locations to values at cone locations
%
% Description
%   Implements the displacement of values from RGC to cone locations.
%
% Inputs:
%   valuesAtRGCSoma       - A 1 x p vector of values at each of the p
%                           eccentricity locations.
%   displacementDegVisual - A 1 x p vector, where p is the number of
%                           eccentricity positions modeled. The values are
%                           the displacement in visual degrees of the RGC
%                           soma away from the fovea.
%   supportPosDegVisual   - A 1 x p vector that contains the locations in
%                           eccentricity in visual degrees for each of
%                           the valuesAtRGCSoma and displacementDegVisual
%
% Optional key/value pairs:
%  'sampleResolutionDegVisual' - The calculations are performed across
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
p.addRequired('displacementDegVisual',@isnumeric);
p.addRequired('supportPosDegVisual',@isnumeric);

% Optional anaysis params
p.addParameter('sampleResolutionDegVisual',0.01,@isnumeric);

% parse
p.parse(valuesAtRGCSoma, displacementDegVisual, supportPosDegVisual, varargin{:})

% Implement the transformation

% We first identfy the destination (in eccentricities) at each of the cone
% locations for each of the RGC locations (rounded to the closest support
% position)
destination=round(supportPosDegVisual-displacementDegVisual,ceil(-log10(p.Results.sampleResolutionDegVisual)));
uniqueDestinations = unique(destination);
destinationIdx=arrayfun(@(x) find((abs(supportPosDegVisual-x))<1e-6,1), uniqueDestinations, 'UniformOutput', false);
validDestinationIndices = cellfun(@(x) ~isempty(x), destinationIdx);
valuesAtConesAtValidDestinations = arrayfun(@(x) sum(valuesAtRGCSoma(destination==x)), uniqueDestinations(validDestinationIndices));
valuesAtCones=zeros(size(valuesAtRGCSoma));
valuesAtCones(unique(cell2mat(destinationIdx(validDestinationIndices))))=valuesAtConesAtValidDestinations;

end