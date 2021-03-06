function [ mRGCDensitySqDegVisual, midgetFraction ] = transformRGCToMidgetRGCDensityDrasdo( regularSupportPosDegVisual, rgcDensitySqDegVisual, varargin )
% Transform a vector of RGC densities to midget RGC densities (Drasdo)
%
% Description:
%   This function transforms a vector of RGC density measurements into
%   corresponding midget RGC density measurements. The eccentricity
%   location of the measurements are not used explicitly in the
%   calculation, instead the cumulative sum of the RGC density.
%
%   This transformation is under the control of four parameters. The first
%   is the f0 value (from Eq 7 of Watson 2014 JoV) that defines the midget
%   fraction at the fovea. The next three are parameters of a reciprocal
%   function.
%
%   If x is the proportion of the cumulative RGC density function at a
%   given location, then:
%
%       midgetFraction =  f0 - [(1./(a+(b.* log10(x) )))+c]
%
%   where f0 is the maximum midget fraction found at the fovea, given by
%   Watson / Drasdo as 0.8928.
%
% Inputs:
%   regularSupportPosDegVisual - A 1 x p vector that contains the
%                           eccentricity in visual degrees at which the
%                           model was evaluated along each meridian
%   rgcDensitySqDegVisual - A 1 x p vector of RGC density at each
%                           eccentricity location in units of cells /
%                           degree^2
%
% Optional key/value pairs:
%  'referenceEccenDegVisual' - The reference eccentricity for the
%                           proportion of the cumulative RGC density. The
%                           proportion function will have a value of unity
%                           at this point.
%  'watsonEq8_f0'         - The midget fraction assigned to the fovea
%  'linkingFuncParams'    - Parameters of the reciprocal fit the defines
%                           the transformation.
%
% Outputs:
%   mRGCDensitySqDegVisual - Midget RGC density at each eccentricity
%                           location
%   midgetFraction        - The midget fraction at each eccentricity
%                           location
%


%% Parse input and define variables
p = inputParser;

% required input
p.addRequired('regularSupportPosDegVisual',@isnumeric);
p.addRequired('rgcDensitySqDegVisual',@isnumeric);

% Optional anaysis params
p.addParameter('totalRetinalRGCs',1.5e6,@isscalar);
p.addParameter('watsonEq8_f0',0.8928,@isnumeric);
p.addParameter('linkingFuncParams',[5.3136   -9.9397   -0.0111],@isnumeric);

% parse
p.parse(regularSupportPosDegVisual, rgcDensitySqDegVisual, varargin{:})


%% House keeping and setup

% Define a three-parameter reciprocal function that will be used to fit the
% modeled relationship
recipFunc = fittype('(1./(a+(b.*x)))+c','independent','x','dependent','y');

% Obtain the cumulative RGC function
RGC_ringcount = calcRingCumulative(regularSupportPosDegVisual,rgcDensitySqDegVisual);

% Calculate a proportion of the cumulative RGC density counts, relative
% to the total number of RGCs in the retina
propRGC_ringcount=RGC_ringcount./p.Results.totalRetinalRGCs;

% Because we are going to be working with a log transform, set any zero
% proportion values to the minimum, non-zero proportion value
zeroPoints=find(propRGC_ringcount==0);
if ~isempty(zeroPoints)
    propRGC_ringcount(zeroPoints)=min(propRGC_ringcount(find(propRGC_ringcount~=0)));
end

% Calculate the midgetFraction based upon the propRGC_ringcount
midgetFraction = p.Results.watsonEq8_f0-recipFunc(p.Results.linkingFuncParams(1),p.Results.linkingFuncParams(2),p.Results.linkingFuncParams(3),log10(propRGC_ringcount));

% Scale the rgcDensity by the midget fraction
mRGCDensitySqDegVisual = rgcDensitySqDegVisual .* midgetFraction;

% NaNs can happen; set them to zero
badIdx = find(isnan(mRGCDensitySqDegVisual));
if ~isempty(badIdx)
    mRGCDensitySqDegVisual(badIdx)=0;
end

end % function


