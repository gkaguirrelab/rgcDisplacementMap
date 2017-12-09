function cumulativeCountsPerRing = calcRingCumulative(regularSupportPosDegRetina, densityFunction)
% Calculate the cumulative number of cells within increasing areal rings
%
% Description:
%   For each of the positions defined in regularSupportPosDegRetina, we
%   calculate the area of a ring that is centered on the fovea (0 degrees
%   retinal eccentricity) and has a width equal to the distance between
%   the support positions. Across these rings of expanding area, we
%   calculate the cumulative number of cells at each position. This vector
%   is used subsequently to identify points on the retina that have matched
%   cumulative sums so as to derive retinal displacement.
%
%   If we wished to know the actual number of cells in each "sector" of
%   retina along a meridian, then the cumulativeCountsPerRing could be
%   defined by the number of polarAngle divisions that are being examined
%   in the model. In practice, we don't bother, as this multiplicative
%   correction would be applied equally to the mRF and mRGC ring cumulative
%   measures and thus would not alter the model output.
%
% Inputs:
%   regularSupportPosDegRetina - A 1 x p vector that contains the
%                           eccentricity in retinal degrees at which the
%                           model was evaluated along each meridian
%   densityFunction       - A 1 x p vector containing the cell (or
%                           receptive field) density at each of the
%                           locations defined in the support.
%
% Outputs:
%   cumulativeCountsPerRing - A 1 x p vector containing the cumulative cell
%                           (or RF) counts per areal ring at each of the
%                           regular support positions
%

% We extend the locatio support by one element to allow us to conduct a
% diff operation without shrinking our array.
extendedRegularSupportPosDegRetina = ...
    [regularSupportPosDegRetina (regularSupportPosDegRetina(end)+diff(regularSupportPosDegRetina(1:2)))];

% We obtain the area of each ring that is centered on the fovea
% (eccentricity 0) and has a width equal to the spacing between the regular
% support locations
ringArea = diff(extendedRegularSupportPosDegRetina.^2 * pi);

% The number of cells / RFs in each ring is given by the product of the
% density function and the ring area at each support position. We then
% calculate and return the cumulative sum of this vector.
cumulativeCountsPerRing = cumsum(densityFunction.*ringArea);

end