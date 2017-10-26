function [ mRFDensitySqDeg, mRFtoConeDensityRatio ] = transformConeToMidgetRFDensity( coneDensitySqDeg, varargin )
% transformConeToMidgetRFDensity( coneDensitySqDeg, varargin )
%
% This function transforms a vector of cone density measurements into
% corresponding midget receptive field density measurements. The
% eccentricty or polar angle location of the measurement is not used
% explicitly in the calculation.
%
% If x is the proportion of the cone density relative to the maxium cone
% density at the fovea, then the ratio of midget receptive fields to cones
% is given by:
%
%	mRFtoConeDensityRatio = minMidgetRGCToConeRatio+(maxMidgetRGCToConeRatio-minMidgetRGCToConeRatio)./(1+(x./inflect).^slope)
%
% INPUT
%   coneDensitySqDeg - a vector of cone densities
%
% OUTPUT
%   mRFDensitySqDeg - a vector of midget receptive field density values
%   mRFtoConeDensityRatio - a vector of midgetRF to cone ratio values
%
% OPTIONS
%   maxConeDensity - The maximum cone density at the fovea (counts / deg^2).
%       The default value is derived from Curcio 1990. If set to empty, the
%       maximum value from coneDensitySqDeg is used.
%   minMidgetRGCToConeRatio - The minimum value of the mRF:cone density ratio. Set to zero
%       as the functions appear to asymptote close to this value.
%   maxMidgetRGCToConeRatio - The maximuim value of the mRF:cone density ratio.
%   linkingFuncParams - parameters for the logisitic fit, corresponding to the
%       slope and inflection point. The default values are those found by
%       fitting the Curcio data from the four meridiansl

%% Parse input and define variables
p = inputParser;

% required input
p.addRequired('coneDensitySqDeg',@isnumeric);

% Optional anaysis params
p.addParameter('maxConeDensity',1.4806e+04,@(x)(isempty(x) | isnumeric(x)));
p.addParameter('minMidgetRGCToConeRatio',-0.1,@isnumeric);
p.addParameter('maxMidgetRGCToConeRatio',2,@isnumeric);
p.addParameter('linkingFuncParams',[5.9861, -1.0636],@isnumeric);

% parse
p.parse(coneDensitySqDeg, varargin{:})

% Set the maxConeDensity value
if isempty(p.Results.maxConeDensity)
    maxConeDensity = max(coneDensitySqDeg);
else
    maxConeDensity = p.Results.maxConeDensity;
end

% Define a four-parameter logistic function that will be used to fit the
% modeled relationship. Two of the parameters (max and min asymptote) are
% locked by the passed parameter
logisticFunc = fittype( @(slope,inflect,minMidgetRGCToConeRatio,maxMidgetRGCToConeRatio,x) minMidgetRGCToConeRatio+(maxMidgetRGCToConeRatio-minMidgetRGCToConeRatio)./(1+sign(x./inflect).*abs((x./inflect).^slope)), ...
    'independent','x','dependent','y','problem',{'minMidgetRGCToConeRatio','maxMidgetRGCToConeRatio'});

% Define the x-axis as the log10 of the proportion of max cone density
x = log10(coneDensitySqDeg ./ maxConeDensity)';

% Obtain the midgetRF : cone ratio from the logisitc function
mRFtoConeDensityRatio = ...
    logisticFunc(p.Results.linkingFuncParams(1), p.Results.linkingFuncParams(2), p.Results.minMidgetRGCToConeRatio, p.Results.maxMidgetRGCToConeRatio, x);

% Calculate the mRF density
mRFDensitySqDeg = coneDensitySqDeg .* mRFtoConeDensityRatio';


end % function


