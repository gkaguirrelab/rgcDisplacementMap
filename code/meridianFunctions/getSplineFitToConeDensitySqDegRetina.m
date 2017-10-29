function [fitConeDensitySqDegRetina] = getSplineFitToConeDensitySqDegRetina(polarAngle, varargin)
% getSplineFitToConeDensitySqDegRetina(angle)
%
% This routine returns a fit to cone density data. Raw cone density data
% are taken from Curcio et al (1990) unless otherwise specified.
% Fits to the cardinal meridians are obtained. Interpolation over the
% parameters of the fit are used to produce a function that returns cone
% density for an arbitrary meridian angle.
%
% Inputs:
%   polarAngle - The desired angle of the density function on the retinal
%                 field. (0=nasal;90=superior;180=temporal;270=inferior)
% Outputs:
%   fitConeDensitySqDegRetina - handle of a fitting function that returns
%       cone density values for the specified meridian angle across retinal
%       eccentricity in polarAngle
%
% Options:

%% Parse input and define variables
p = inputParser;

% required input
p.addRequired('polarAngle',@isnumeric);

% Optional analysis params
p.addParameter('cardinalMeridianAngles',[0, 90, 180, 270],@isnumeric);
p.addParameter('splineKnots',15,@isnumeric);
p.addParameter('splineOrder',4,@isnumeric);
p.addParameter('coneDensityDataFileName', [], @(x)(isempty(x) | ischar(x)));

% parse
p.parse(polarAngle,varargin{:})

%% sanity check the input
if p.Results.cardinalMeridianAngles ~= 0
    error('This routine assumes that the first meridian is polar angle = 0');
end

% Loop across the cardinal meridians and combine the cone data. Use this
% aggregate data to define the knots for the spline fit.
aggregatePosition=[];
aggregateDensity=[];
for mm=1:length(p.Results.cardinalMeridianAngles)
    % get the raw density measurements
    if isempty(p.Results.coneDensityDataFileName)
        % load the empirical cone density measured by Curcio
        [coneDensitySqDegRetina, coneNativeSupportPosDegRetina] = ...
            loadRawConeDensityByEccen(p.Results.cardinalMeridianAngles(mm));
    else
        % or the passed cone density measurement
        [coneDensitySqDegRetina, coneNativeSupportPosDegRetina] = ...
            loadRawConeDensityByEccen(p.Results.cardinalMeridianAngles(mm), ...
            'coneDensityDataFileName', p.Results.coneDensityDataFileName);
    end
    
    % remove nan values
    isvalididx=find(~isnan(coneDensitySqDegRetina));
    coneNativeSupportPosDegRetina = coneNativeSupportPosDegRetina(isvalididx);
    coneDensitySqDegRetina = coneDensitySqDegRetina(isvalididx);
    aggregatePosition = [aggregatePosition coneNativeSupportPosDegRetina];
    aggregateDensity = [aggregateDensity coneDensitySqDegRetina];
end

% perform a least-squares spline fit with the specified knots and order
BformSplineFit=spap2(p.Results.splineKnots, p.Results.splineOrder, aggregatePosition', aggregateDensity');
knots = BformSplineFit.knots;

% Loop across the cardinal meridians again, and now perform the spline fit
% with the specified knots
for mm=1:length(p.Results.cardinalMeridianAngles)
    % load the empirical cone density measured by Curcio
    [coneDensitySqDegRetina, coneNativeSupportPosDegRetina] = loadRawConeDensityByEccen(p.Results.cardinalMeridianAngles(mm));
    % remove nan values
    isvalididx=find(~isnan(coneDensitySqDegRetina));
    coneNativeSupportPosDegRetina = coneNativeSupportPosDegRetina(isvalididx);
    coneDensitySqDegRetina = coneDensitySqDegRetina(isvalididx);
    % Perform the spline fit
    BformSplineFit=spap2(knots, p.Results.splineOrder, coneNativeSupportPosDegRetina', coneDensitySqDegRetina');
    % convert from B-form to piecewise polynomial form
    ppFormSplineFits{mm} = fn2fm(BformSplineFit,'pp');
end

% We will now assemble a new coefficient matrix for the spline fit that is
% an interpolation between the coefficient matrices obtained for the
% cardinal meridians. This is done by going through each element of the
% coefficient matrix and obtaining the four values, one from each meridian.
% The value from polar angle 0 is replicated and assigned a polar angle of
% 360; this allows the interpolation to wrap-around. We then fit a 'pchip'
% (Piecewise Cubic Hermite Interpolating Polynomial) to the set of 5 values
% and interpolate to the called-for, intermediate polar angle.

% make an empty coefficient matrix
interpCoefs = zeros(size(ppFormSplineFits{1}.coefs));

% replicate the values for polar angle zero at polar angle 360 to allow a
% wrap-around fit
angleBase=[p.Results.cardinalMeridianAngles 360];

% loop over each element of the coefficient matrix
for cc=1:numel(interpCoefs)
    % obtain the set of coefficients from across the meridians
    thisCoefByMeridian = cellfun(@(x) x.coefs(cc), ppFormSplineFits);
    % duplicate the first coefficient (assuming it is polar angle zero)
    thisCoefByMeridian = [thisCoefByMeridian thisCoefByMeridian(1)];
    % fit the pchip model
    coefByAngleFit = fit(angleBase',thisCoefByMeridian','pchipinterp');
    % evaluate the fitted model at the called for polar angle and place the
    % returned, interpolated coefficient value into the matrix we are
    % building
    interpCoefs(cc) = coefByAngleFit(polarAngle);
end

% Assemble our interpolated piecewise polynomial fit structure
ppFormSplineInterp = ppFormSplineFits{1};
ppFormSplineInterp.coefs = interpCoefs;

% Create an anonymous function using the interpolated spline. This function
% will return cone density as a function of eccentricity in degrees
fitConeDensitySqDegRetina = @(supportPosDeg) fnval(ppFormSplineInterp,supportPosDeg');

end % function

