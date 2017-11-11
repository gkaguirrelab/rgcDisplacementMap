function valuesAtCones = applyDisplacement(valuesAtRGCSoma, displacementDegRetina, supportPosDegRetina, varargin)
% applyDisplacement
%
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

destination=round(supportPosDegRetina-displacementDegRetina,ceil(-log10(p.Results.sampleResolutionDegreesRetina)));
uniqueDestinations = unique(destination);
destinationIdx=arrayfun(@(x) find((abs(supportPosDegRetina-x))<1e-6,1), uniqueDestinations, 'UniformOutput', false);
validDestinationIndices = cellfun(@(x) ~isempty(x), destinationIdx);
valuesAtConesAtValidDestinations = arrayfun(@(x) sum(valuesAtRGCSoma(destination==x)), uniqueDestinations(validDestinationIndices));
valuesAtCones=zeros(size(valuesAtRGCSoma));
valuesAtCones(unique(cell2mat(destinationIdx(validDestinationIndices))))=valuesAtConesAtValidDestinations;

end