function [ mRGCDensitySqDeg, midgetFraction ] = transformRGCToMidgetRGCDensity( regularSupportPosDeg, rgcDensitySqDeg, varargin )
% transformRGCToMidgetRGCDensity( regularSupportPosDeg, rgcDensitySqDeg, varargin )
%
% This function transforms a vector of RGC density measurements into
% corresponding midget RGC density measurements. The eccentricity location
% of the measurements are not used explicitly in the calculation, instead
% the cumulative sum of the RGC density.
%
% This transformation is under the control of four parameters. Two of these
% are the max and min midget ratios found in the fovea and periphery,
% respectively; these are treated as fixed parameters of the model. The
% other two parameters are the slope and inflection point of the logisitic
% fit.
%
% If x is the proportion of the cumulative RGC density function at a given
% location, then:
%
%   midgetFraction = minMidgetFractionRatio+(maxMidgetFractionRatio-minMidgetFractionRatio)./(1+sign(x./inflect).*abs((x./inflect).^slope))
%
% where inflect and slope are free paramters of the logisitic fit, and
% maxMidgetFractionRatio and minMidgetFractionRatio are the values of the midget fraction at the fovea
% and periphery, respectivey.
%
% INPUT
%   regularSupportPosDeg - The eccentricity values that support
%       rgcDensitySqDeg.
%   rgcDensitySqDeg - RGC density at each eccentricity location in units of
%      cells / degree^2
%
% OUTPUT
%   mRGCDensitySqDeg - midget RGC density at each eccentricity location
%   midgetFraction - the midget fraction at each eccentricity location
%
% OPTIONS
%   referenceEccen - the reference eccentricity for the proportion of
%       the cumulative RGC density. The proportion function will have a
%       value of unity at this point.
%   maxMidgetFractionRatio - The midget fraction assigned to the fovea
%   minMidgetFractionRatio - The midget fraction assigned to the far periphery
%   linkingFuncParams - the slope and inflection parameters that are used in
%       the logisitic function.
%   verbose - Controls text output to console
%   makePlots - Do we make a figure?

%% Parse input and define variables
p = inputParser;

% required input
p.addRequired('regularSupportPosDeg',@isnumeric);
p.addRequired('rgcDensitySqDeg',@isnumeric);

% Optional anaysis params
p.addParameter('referenceEccen',15,@isnumeric);
p.addParameter('minMidgetFractionRatio',0.41,@isnumeric);
p.addParameter('maxMidgetFractionRatio',0.85,@isnumeric);
p.addParameter('linkingFuncParams',[5 1],@isnumeric);

% parse
p.parse(regularSupportPosDeg, rgcDensitySqDeg, varargin{:})


%% Perform the calculations

% Define a four-parameter logistic function that will be used to fit the
% modeled relationship. Two of the parameters (max and min asymptote) are
% locked by the passed parameter
logisticFunc = fittype( @(slope,inflect,minMidgetFractionRatio,maxMidgetFractionRatio,x) minMidgetFractionRatio+(maxMidgetFractionRatio-minMidgetFractionRatio)./(1+sign(x./inflect).*abs((x./inflect).^slope)), ...
    'independent','x','dependent','y','problem',{'minMidgetFractionRatio','maxMidgetFractionRatio'});

% Obtain the cumulative RGC function
RGC_ringcount = calcCumulative(regularSupportPosDeg,rgcDensitySqDeg);

% Find the index position in the regularSupportPosDeg that is as close
% as possible to the referenceEccen
[ ~, refPointIdx ] = min(abs(regularSupportPosDeg-p.Results.referenceEccen));

% Calculate a proportion of the cumulative RGC density counts, relative
% to the reference point (which is assigned a value of unity)
propRGC_ringcount=RGC_ringcount./RGC_ringcount(refPointIdx);

% Calculate the midgetFraction based upon the propRGC_ringcount
midgetFraction = logisticFunc(p.Results.linkingFuncParams(1), p.Results.linkingFuncParams(2), p.Results.minMidgetFractionRatio, p.Results.maxMidgetFractionRatio, propRGC_ringcount);

% Scale the rgcDensity by the midget fraction
mRGCDensitySqDeg = rgcDensitySqDeg .* midgetFraction;

% NaNs can happen; set them to zero
badIdx = find(isnan(mRGCDensitySqDeg));
if ~isempty(badIdx)
    mRGCDensitySqDeg(badIdx)=0;
end

end % function
