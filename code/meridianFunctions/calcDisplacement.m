function displaceInDegRetina = calcDisplacement(regularSupportPosDegRetina, mRF_RingCumulative, mRGC_RingCumulative)
% Calculate RGC displacement using RF and RGC cumulative functions
%
% Description:
%   This function identifies points along the regularSupportPosDegRetina at
%   which the midget RGC and midget receptive field (RF) cumulative
%   functions have matching values. The distance (in retinal degrees)
%   between these matching points is calculated and returned.
%
% Inputs:
%   regularSupportPosDegRetina - A 1 x p vector that contains the
%                           eccentricity in retinal degrees at which the
%                           model was evaluated along each meridian
%   mRF_RingCumulative    - A 1 x p vector that contains the cumulative
%                           number of RFs at each regularSupportPos-
%                           DegRetina position within an areal "ring".
%   mRGC_RingCumulative   - A 1 x p vector that contains the cumulative
%                           number of RGC cells at each regularSupportPos-
%                           DegRetina position within an areal "ring".
%
% Outputs:
%   displaceInDegRetina   - A 1 x p vector that contains, for each of the
%                           retinal positions defined in the support,
%                           the distance (in retinal degrees) that RGC cell
%                           bodies at that location would need to be
%                           displaced towards the retina to place them at
%                           the same location as their receptive fields.
%


% Determine the sample resolution by a difference operation
tmp = diff(regularSupportPosDegRetina);
sampleResolutionDegreesRetina = tmp(1);

% Measure the displacement (in degrees). First, for each cumulative RGC
% density value, identify the array index of the last RF density value
% that is less than or equal to the RGC value.
displaceInSamples=arrayfun(@(x) find(mRF_RingCumulative<=x,1,'last')+1, mRGC_RingCumulative,'UniformOutput',false);

% Handle the case of the initial RF values, all of which are larger than
% the first few RGC values
displaceInSamples(find(mRF_RingCumulative(1)>mRGC_RingCumulative))={1};

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