function [ fitParams, figHandle ] = developMidgetRGCFractionModel( varargin )
% developMidgetRGCFractionModel( varargin )
%
% Our goal is to be able to derive the fraction of RGCs that are midget
% RGCs at a point in the retina, but without reference to the polar angle
% or eccentricity of the measurement.
%
% To solve this problem, we first take the cumulative of the RGC density
% function along a meridian. We then express this function as a proportion
% of the cumuluative density value at any point to the cumulative density
% value at a specific eccentricity location (e.g., at 30 degrees). Thus,
% the proportional cumulative density of RGCs is equal to unity (1) at the
% reference eccentricity. The value of the reference eccentricity should be
% sufficiently large so that the RGC density function is in the
% monotonically decreasing phase.
%
% Next, the log10 of the proportional, cumulative RGC density function is
% related to the predicted midget fraction for that location, as derived by
% eq 7 of Watson 2014 JOV (which itself is taken from Drasdo 2007). We
% observe that a three-component reciprocal function fits this relationship
% well.
%
% We note here a conceptual inconsistency in our approach. Our goal is to
% model the proportion of RGCs that are midget RGCs at the retinal location
% that contains these cell bodies. We use the Watson equation (via Drasdo)
% to set the "ground truth" of the midget fraction, but this function
% expresses the midget fraction at the receptive field locations of these
% cells. We would ideally have a measurement of the fraction of midget RGC
% soma at each of many eccentricities in the human retina, but this is not
% available.
%
% We obtain the fit parameters for the midget fraction expression for each
% of the four meridians for which we have data, and take the median fit
% parameters across merdians. This gives us a functional form that can be
% used to convert RGC density to midget RGC density knowing only the
% proportion of the cumulative RGC density function at the location for
% which the measurement was made. If we term this proportion value x, then
% the relationship is:
%
%   midgetFraction =  f0 - [(1./(a+(b.* log10(x) )))+c]
%
% where f0 is the maximum midget fraction found at the fovea, given by
% Watson / Drasdo as 0.8928.
%
% OUTPUT
%   fitParams - The median values of [a b c] across the four merdians
%   figHandle - Handle to a figure showing the fits. Empty if no plotting
%      requested.
%
% OPTIONS
%   referenceEccen - the reference eccentricity for the proportion of
%       the cumulative RGC density. The proportion function will have a
%       value of unity at this point. We use 15° here for the practical
%       reason that this is the maximum extent for which we have OCT
%       measurements of the RGC layer thickness, and we wish in the future
%       to model such data using these functions.
%   supportResolutionDegrees - the resolution (in degrees) for which the
%       calculations are made.
%   supportEccenMaxDegrees - the maximum eccentricity used for modeling
%   meridianNames - Cell array of the text string names of the meridia
%   meridianAngles - Polar angle values assigned to the meridians
%   watsonEq8_f0 - The midget fraction assigned to the fovea
%   watsonEq8_rm - decay parameter controlling the midget fraction function
%   recipFitStartPoint - initial values used for the fit of the
%       three-parameter reciprocal function. It is critical that the second
%       parameter (b) has a negative start point to fit this inverted
%       reciprocal.
%   makePlots - Do we make a figure?

%% Parse input and define variables
p = inputParser;

% Optional anaysis params
p.addParameter('referenceEccen',15,@isnumeric);
p.addParameter('supportResolutionDegrees',0.01,@isnumeric);
p.addParameter('supportEccenMaxDegrees',30,@isnumeric);
p.addParameter('meridianNames',{'Nasal' 'Superior' 'Temporal' 'Inferior'},@iscell);
p.addParameter('meridianAngles',[0, 90, 180, 270],@isnumeric);
p.addParameter('meridianSymbols',{'.','x','o','^'},@cell);
p.addParameter('minRatio',0.45,@isnumeric);
p.addParameter('maxRatio',0.95,@isnumeric);
p.addParameter('logitFitStartPoint',[-1 1],@isnumeric);

% Optional display params
p.addParameter('makePlots',true,@islogical);

% parse
p.parse(varargin{:})


%% House keeping and setup

% Define a four-parameter logistic function that will be used to fit the
% modeled relationship. Two of the parameters (max and min asymptote) are
% locked by the passed parameter
logisticFunc = fittype( @(slope,inflect,minRatio,maxRatio,x) minRatio+(maxRatio-minRatio)./(1+sign(x./inflect).*abs((x./inflect).^slope)), ...
    'independent','x','dependent','y','problem',{'minRatio','maxRatio'});

% Define the regular-spaced eccentricity support over which we will model
% the anatomical retinal functions
regularSupportPosDeg = ...
    0:p.Results.supportResolutionDegrees:p.Results.supportEccenMaxDegrees;

% Prepare a figure if requested
if p.Results.makePlots
    figHandle=figure;
else
    figHandle=[];
end

%% Loop over the meridians

for mm = 1:length(p.Results.meridianAngles)
    
    % Load the RGC Density Data from Curcio and Allen 1990
    [ RGCDensitySqDeg, nativeSupportPosDeg ] = getCurcioRGCDensityByEccen( p.Results.meridianAngles(mm) );
    
    % remove nan values
    isvalididx=find(~isnan(RGCDensitySqDeg)  );
    nativeSupportPosDeg = nativeSupportPosDeg(isvalididx);
    RGCDensitySqDeg = RGCDensitySqDeg(isvalididx);
    
    % Fit a spline to the RGC density data
    RGCDensityFit = fit(nativeSupportPosDeg',RGCDensitySqDeg','smoothingspline', 'Exclude',find(isnan(RGCDensitySqDeg)),'SmoothingParam', 1);
    
    % Obtain the cumulative RGC function
    RGC_ringcount = calcCumulative(regularSupportPosDeg,RGCDensityFit(regularSupportPosDeg)');
    
    % Find the index position in the regularSupportPosDeg that is as close
    % as possible to the referenceEccen
    refPointIdx= ...
        find((regularSupportPosDeg-p.Results.referenceEccen)== ...
        min(abs(regularSupportPosDeg-p.Results.referenceEccen)));
    % Calculate a proportion of the cumulative RGC density counts, relative
    % to the reference point (which is assigned a value of unity)
    propRGC_ringcount=RGC_ringcount./RGC_ringcount(refPointIdx);
        
    % Obtain the Dacey midget fraction as a function of eccentricity
    midgetFractionByEccen = calcDaceyMidgetFractionByEccen(regularSupportPosDeg)';
    
    % Fit the logisitc model that relates proportionRGC to midget
    % fraction. First, define a weight function to lock the max value
    weights=ones(1,length(propRGC_ringcount));
    weights(1)=1e5;
    
    % Perform the fit and save the param values.
    % The x values proportionRGC
    % The y values are the midget fraction
    
    % turn off some warnings that are produced by the fit and plot
    warningState = warning;
    warning('off','curvefit:fit:complexXusingOnlyReal');
    warning('off','MATLAB:plot:IgnoreImaginaryXYPart');
    
        % Perform the logistic fit. Note that the max and min asymptote are
    % pinned by the passed parameters
    logitFit = fit(propRGC_ringcount',midgetFractionByEccen',logisticFunc, ...
        'problem',{p.Results.minRatio, p.Results.maxRatio}, ...
        'StartPoint',p.Results.logitFitStartPoint, ...
        'Lower',[-10,-10],'Upper',[10,10] );
    
    % Distribute the logitFit params to a loop variable
    fitParams(mm,1)=logitFit.slope;
    fitParams(mm,2)=logitFit.inflect;
    
    % Add the data to the figure
    if p.Results.makePlots
        plot(propRGC_ringcount,midgetFractionByEccen,p.Results.meridianSymbols{mm},'color',[.8 .8 .8])
        hold on
        ylim([0.4 1]);
        xlim([0 2]);
        xlabel('proportion cumulative RGC density count at reference eccen');
        ylabel('midget fraction');
    end % if we are plotting
    
    % restore the saved warning state
    warning(warningState);
    
end % loop over meridians

% Now calculate the median parm values across the meridians and add a fit
% line to the plot
fitParams=median(fitParams);

if p.Results.makePlots
    xFit=0:.1:2;
    plot( xFit, ...
        logisticFunc(fitParams(1), fitParams(2), p.Results.minRatio, p.Results.maxRatio, xFit),'-r')
    legend({p.Results.meridianNames{:} 'fit'},'Location','southwest');
    title('midget fraction as a function of relative RGC cumulative density');
    drawnow
end


end % function


