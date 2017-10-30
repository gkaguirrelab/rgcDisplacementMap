function displaceInDegRetina = calcDisplacement(regularSupportPosDegRetina, countPerRingRGC, countPerRingRF)

% Determine the sample resolution by a difference operation
tmp = diff(regularSupportPosDegRetina);
sampleResolutionDegreesRetina = tmp(1);

% Measure the displacement (in degrees). First, for each cumulative RGC
% density value, identify the array index of the last RF density value
% that is less than or equal to the RGC value.
displaceInSamples=arrayfun(@(x) find(countPerRingRF<=x,1,'last')+1, countPerRingRGC,'UniformOutput',false);

% Handle the case of the initial RF values, all of which are larger than
% the first few RGC values
displaceInSamples(find(countPerRingRF(1)>countPerRingRGC))={1};

% Now some array operations to get these values out of cells and in to a
% numeric vector
emptyCells = find(cellfun(@(x) isempty(x), displaceInSamples));
displaceInSamples(emptyCells)={NaN};
displaceInSamples=cell2mat(displaceInSamples(cellfun(@(x) ~isempty(x), displaceInSamples)));

% The displacement in array samples is the array index of each cumulative
% RGC density measurement, minus the array index of the matching RF density
% measurement. This difference is then multiplied by the sample resolution
% to get the displacement in degrees.
displaceInDegRetina = ((1:1:length(displaceInSamples))-displaceInSamples ) * sampleResolutionDegreesRetina;

% Zero out negative values after the point of convergence
displaceInDegRetina(find(displaceInDegRetina < 0,1):end)=0;

end