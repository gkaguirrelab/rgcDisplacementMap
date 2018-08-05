function [fitRGCDensitySqDegVisual, figHandle] = getSplineFitToRGCDensitySqDegVisual(polarAngle, varargin)
% Returns a fit to RGC density data
%
% Syntax:
%  [fitRGCDensitySqDegVisual, figHandle]] = getSplineFitToRGCDensitySqDegVisual(polarAngle)
%
% Description:
%   This routine returns a fit to RGC density data. Raw RGC density data
%   are taken from Curcio et al (1990) unless otherwise specified. Fits to
%   the cardinal meridians are obtained. Interpolation over the parameters
%   of the fit are used to produce a function that returns RGC density for
%   an arbitrary meridian angle.
%
% Inputs:
%   polarAngle            - The desired angle of the density function on 
%                           the retinal field. (0=nasal;90=superior;
%                           180=temporal;270=inferior)
%
% Optional key/value pairs:
%  'cardinalMeridianAngles' - The polar angles corresponding to the
%                           cardinal medians
%  'splineKnots'          - The number of spline knots to use in the fit
%  'splineOrder'          - The polynomial order of the spline
%  'rgcDensityDataFileName' - The filename of the RGC density file to be
%                           passed to loadRawRGCDensityByEccen. If set to
%                           empty, then the default setting in the load
%                           routine will be used.
%  'meridiansForKnotDefinition' - We find that the knots must be defined
%                           using only the 180 and 270 degree meridians to
%                           produce knots that work for all meridians. This
%                           is a hack. The indices set here point to these
%                           poar angles in the cardinalMeridianAngles
%                           vector.
%  'makePlots'            - Do we make a figure?
%
% Outputs:
%	fitRGCDensitySqDegVisual - handle of a fitting function that returns
%                           rgc density values for the specified meridian
%                           angle across retinal eccentricity in polarAngle
%   figHandle             - Handle to a figure showing the fits. Empty if 
%                           no plotting requested.
%
% Examples:
%{
    fitRGCDensitySqDegVisual = getSplineFitToRGCDensitySqDegVisual(45, 'makePlots', true);
%}

%% Parse input and define variables
p = inputParser;

% required input
p.addRequired('polarAngle',@isnumeric);

% Optional analysis params
p.addParameter('cardinalMeridianAngles',[0, 90, 180, 270],@isnumeric);
p.addParameter('cardinalMeridianNames',{'nasal','superior','temporal','inferior'},@iscell);
p.addParameter('cardinalMeridianPlotColors',{'r','b','g','k'},@iscell);
p.addParameter('splineKnots',18,@isnumeric);
p.addParameter('splineOrder',4,@isnumeric);
p.addParameter('rgcDensityDataFileName', [], @(x)(isempty(x) | ischar(x)));
p.addParameter('meridiansForKnotDefinition',[1,2,3,4],@isnumeric);

% Optional display params
p.addParameter('makePlots',false,@islogical);

% parse
p.parse(polarAngle,varargin{:})


%% sanity check the input
if p.Results.cardinalMeridianAngles ~= 0
    error('This routine assumes that the first meridian is polar angle = 0');
end

% Prepare a figure if requested
if p.Results.makePlots
    figHandle=figure;
    subplot(1,2,1);
else
    figHandle=[];
end

% Loop across the cardinal meridians and combine the RGC data. Use this
% aggregate data to define the knots for the spline fit.
aggregatePosition=[];
aggregateDensity=[];
for mm=1:length(p.Results.cardinalMeridianAngles)
    % get the raw density measurements
    if isempty(p.Results.rgcDensityDataFileName)
        % load the empirical RGC density measured by Curcio
        [rgcDensitySqDegVisual, rgcNativeSupportPosDegVisual] = ...
            loadRawRGCDensityByEccen(p.Results.cardinalMeridianAngles(mm));
    else
        % or the passed RGC density measurement
        [rgcDensitySqDegVisual, rgcNativeSupportPosDegVisual] = ...
            loadRawRGCDensityByEccen(p.Results.cardinalMeridianAngles(mm), ...
            'rgcDensityDataFileName', p.Results.rgcDensityDataFileName);
    end
    % remove nan values
    isvalididx=find(~isnan(rgcDensitySqDegVisual));
    rgcNativeSupportPosDegVisual = rgcNativeSupportPosDegVisual(isvalididx);
    rgcDensitySqDegVisual = rgcDensitySqDegVisual(isvalididx);
    % make a plot if requested
    if p.Results.makePlots
        loglog(rgcNativeSupportPosDegVisual,rgcDensitySqDegVisual,'x','Color',p.Results.cardinalMeridianPlotColors{mm});
        hold on
    end
    % If this meridian is one of the ones to be used, aggregate the values
    if sum(p.Results.meridiansForKnotDefinition==mm)
        % Replace zeros with values close to zero
        rgcNativeSupportPosDegVisual(rgcNativeSupportPosDegVisual==0)=1e-10;
        rgcDensitySqDegVisual(rgcDensitySqDegVisual==0)=1e-10;
        aggregatePosition = [aggregatePosition rgcNativeSupportPosDegVisual];
        aggregateDensity = [aggregateDensity rgcDensitySqDegVisual];
    end
end

% perform a least-squares spline fit with the specified knots and order
BformSplineFit=spap2(p.Results.splineKnots, p.Results.splineOrder, log10(aggregatePosition)', log10(aggregateDensity)');
knots = BformSplineFit.knots;

% Loop across the cardinal meridians again, and now perform the spline fit
% with the specified knots
for mm=1:length(p.Results.cardinalMeridianAngles)
    % load the empirical rgc density measured by Curcio
    [rgcDensitySqDegVisual, rgcNativeSupportPosDegVisual] = ...
        loadRawRGCDensityByEccen(p.Results.cardinalMeridianAngles(mm));
    % remove nan values
    isvalididx=find(~isnan(rgcDensitySqDegVisual));
    rgcNativeSupportPosDegVisual = rgcNativeSupportPosDegVisual(isvalididx);
    rgcDensitySqDegVisual = rgcDensitySqDegVisual(isvalididx);
    % Replace zeros with values close to zero
    rgcNativeSupportPosDegVisual(rgcNativeSupportPosDegVisual==0)=1e-10;
    rgcDensitySqDegVisual(rgcDensitySqDegVisual==0)=1e-10;
    % Perform the spline fit
    BformSplineFit = ...
        spap2(knots, p.Results.splineOrder, log10(rgcNativeSupportPosDegVisual)', log10(rgcDensitySqDegVisual)');
    % convert from B-form to piecewise polynomial form
    ppFormSplineFits{mm} = fn2fm(BformSplineFit,'pp');
    % Add a fit line to the plot
    if p.Results.makePlots
        fitSupport=0:0.1:max(rgcNativeSupportPosDegVisual);
        loglog(fitSupport,10.^fnval(ppFormSplineFits{mm},log10(fitSupport)),'-','Color',p.Results.cardinalMeridianPlotColors{mm});
    end
end

%  Clean up the plot
if p.Results.makePlots    
    xlabel('log10 Eccentricity [deg visual]');
    ylabel('log10 RGC density [counts / deg visual^2]');
    ylim([1 1e4]);
    xlim([1e-1 1e2]);
    legend({p.Results.cardinalMeridianNames{:} 'fit'},'Location','southwest');
    drawnow
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
% will return rgc density as a function of eccentricity in degrees. The
% transpose operations are needed so that that the function returns a row
% vector of density in response to a row vector of eccentricity support.
fitRGCDensitySqDegVisual = @(supportPosDeg) 10.^fnval(ppFormSplineInterp,log10(supportPosDeg)')';

% Show the interpolated meridian
if p.Results.makePlots
    subplot(1,2,2)
    fitSupport=0:0.1:max(rgcNativeSupportPosDegVisual);
    loglog(fitSupport,fitRGCDensitySqDegVisual(fitSupport),'-.k');
    xlabel('log10 Eccentricity [deg visual]');
    ylabel('log10 RGC density [counts / deg visual^2]');
    ylim([1 1e4]);
    xlim([1e-1 1e2]);
end

    

end % function

