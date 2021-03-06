function [ fitParams, figHandle ] = developMidgetRFFractionModel( varargin )
% Find parameters that transform cone density to midget RF density
%
% Syntax:
%  [ fitParams, figHandle ] = developMidgetRFFractionModel()
%
% Description:
%   Examine the relationship between cone density and midget receptive
%   field density. Note that the mRF density reflects both the ratio of
%   convergence of cones onto retinal ganglion cells, as well as the
%   fraction of retinal ganglion cells at any one point that are midgets.
%
%   This function and its companion (developMidgetRGCFractionModel) both
%   incorporate the effect of only a fraction of all RGCs being midget
%   RGCs. However, this function considers that fraction in relation to
%   receptive fields at cone locations, while the companion function
%   expresses the midget fraction in relation to the RGC cell bodies at
%   their displaced positions.
%
%   Watson (2014), citing Kolb & Dekorver (1991), asserts that the ratio of
%   mRFs cones at the fovea is 2:1. That is, each foveal cone provides
%   input to two mRGCs (an "on" and "off" cell). We adopt this constraint
%   as well, although we note that Watson's functions for mRF density have
%   a discontinuity at the fovea, in that his formula has the mRF:cone
%   ratio asymptote at ~1.9, and then jumps to a value of 2 at the fovea.
%
%   For each meridian, We first load the Curcio cone density measurements.
%   These measurements are expressed at positions in visual degrees. At
%   each of these retinal eccentricity values, we also obtain the midgetRF
%   density. This is derived from Eq 8 of Watson 2014.
%
%   We observe that the ratio of midgetRF density to cone density as a
%   function of log10 cone density is sigmoidal. To support generalization
%   across datasets, we express the x-axis as log10 of the proportion of
%   maximum cone density observed in the fovea. We then fit this
%   relationship with the function:
%
%	mRFtoConeDensityRatio = 
%       minMidgetRGCToConeRatio + 
%       (maxMidgetRGCToConeRatio - minMidgetRGCToConeRatio) ./ 
%       (1+(x./inflect).^slope)
%
%   where maxMidgetRGCToConeRatio and minMidgetRGCToConeRatio are locked
%   parameters.
%
%   We adopt a slightly negative minMidgetRGCToConeRatio to provide a
%   better fit of the sigmoid at the far periphery, although we observe
%   that this value is not physiologically meaningful. Effectively, our
%   model returns impossible values when the cone density reaches less than
%   1% of its peak value at the fovea.
%
%   We calculate the fit of this function across each of the four
%   meridians. We observe that our functional form fits the data from the
%   nasal, temporal, and inferior meridians very well. The form for the
%   superior meridian is poorly fit by this model. Therefore, we take the
%   median of the two free fit parameters (slope and inflect) from just the
%   three well-modeled meridians and return these values.
%
% Inputs:
%   none
%
% Optional key/value pairs:
%  'supportEccenMaxDegreesRetina' - The maximum eccentricity used for
%                           modeling. This value should be sufficiently
%                           high so that we are in the asymptote range of
%                           cone density
%  'meridianNames'        - Cell array of the text string names of the
%                           meridians
%  'meridianAngles'       - Polar angle values assigned to the meridians
%  'meridiansIdxToUseForFitParams' - Index of the meridians from which we 
%                           will calculate the median fit param values.
%  'maxConeDensity'       - The maximum cone density at the fovea (counts / 
%                           deg^2). The default value is derived from
%                           Curcio 1990. If set to empty, the maximum value
%                           from coneDensitySqDeg is used.
%  'minMidgetRGCToConeRatio' - The minimum value of the mRF:cone density
%                           ratio.
%  'maxMidgetRGCToConeRatio' - The maximuim value of the mRF:cone density
%                           ratio.
%  'logitFitStartPoint'   - Initial values used for the slope and 
%                           inflection point parameters of the logisic fit.
%                           Hand-crafted after examination of typical data.
%  'makePlots'            - Do we make a figure?
%
% Outputs:
%   fitParams             - The median values of [slope inflect] across the
%                           merdians that were selected to be used
%   figHandle             - Handle to a figure showing the fits. Empty if
%                           no plotting requested.
%
% Examples:
%{
    [ fitParams, figHandle ] = developMidgetRFFractionModel( 'makePlots', true )
%}

%% Parse input and define variables
p = inputParser;

% Optional anaysis params
p.addParameter('supportEccenMaxDegreesRetina',60,@isnumeric);
p.addParameter('meridianNames',{'Nasal' 'Superior' 'Temporal' 'Inferior'},@iscell);
p.addParameter('meridianAngles',[0, 90, 180, 270],@isnumeric);
p.addParameter('meridianSymbols',{'.','x','o','^'},@cell);
p.addParameter('meridiansIdxToUseForFitParams',[1 3 4],@isnumeric);
p.addParameter('maxConeDensitySqDegVisual',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('minMidgetRGCToConeRatio',-.5,@isnumeric);
p.addParameter('maxMidgetRGCToConeRatio',2,@isnumeric);
p.addParameter('logitFitStartPoint',[3,-1],@isnumeric);

% Optional display params
p.addParameter('makePlots',false,@islogical);

% parse
p.parse(varargin{:})


%% House keeping and setup

% Define a four-parameter logistic function that will be used to fit the
% modeled relationship. Two of the parameters (max and min asymptote) are
% locked by the passed parameter
logisticFunc = fittype( @(slope,inflect,minMidgetRGCToConeRatio,maxMidgetRGCToConeRatio,x) minMidgetRGCToConeRatio+(maxMidgetRGCToConeRatio-minMidgetRGCToConeRatio)./(1+sign(x./inflect).*abs((x./inflect).^slope)), ...
    'independent','x','dependent','y','problem',{'minMidgetRGCToConeRatio','maxMidgetRGCToConeRatio'});

% Prepare a figure if requested
if p.Results.makePlots
    figHandle=figure;
else
    figHandle=[];
end

%% Loop over the meridians

for mm = 1:length(p.Results.meridianAngles)
    
    % Load the Cone density Data from Curcio 1990:
    [coneDensitySqDegVisual, nativeSupportPosDegVisual] = loadRawConeDensityByEccen(p.Results.meridianAngles(mm));
    
    % remove nan values
    isvalididx=find(~isnan(coneDensitySqDegVisual));
    nativeSupportPosDegVisual = nativeSupportPosDegVisual(isvalididx);
    coneDensitySqDegVisual = coneDensitySqDegVisual(isvalididx);
    
    % Set the maxConeDensity value
    if isempty(p.Results.maxConeDensitySqDegVisual)
        maxConeDensitySqDegVisual = max(coneDensitySqDegVisual);
    else
        maxConeDensitySqDegVisual = p.Results.maxConeDensitySqDegVisual;
    end
        
    % calculate the mRF density at these eccentricity locations using
    % Watson equation 8.
    [ midgetRFDensitySqDegVisual ] = calcWatsonMidgetRFDensityByEccenDegVisual(nativeSupportPosDegVisual, p.Results.meridianAngles(mm));
        
    % Remove nans and points beyond the modeled eccentricity bound
    isvalididx=find(~isnan(midgetRFDensitySqDegVisual).*~isnan(coneDensitySqDegVisual) .* (nativeSupportPosDegVisual < p.Results.supportEccenMaxDegreesRetina) );
    nativeSupportPosDegVisual = nativeSupportPosDegVisual(isvalididx);
    coneDensitySqDegVisual = coneDensitySqDegVisual(isvalididx);
    midgetRFDensitySqDegVisual = midgetRFDensitySqDegVisual(isvalididx);
    
    % Define the ratio function.
    midgetRFtoConeRatio = (midgetRFDensitySqDegVisual ./ coneDensitySqDegVisual)';
    
    % Define the x-axis as the log10 of the proportion of max cone density
    x = log10(coneDensitySqDegVisual ./ maxConeDensitySqDegVisual)';
        
    % Perform the logistic fit. Note that the max and min asymptote are
    % pinned by the passed parameters
    logitFit = fit(x,midgetRFtoConeRatio,logisticFunc, ...
        'problem',{p.Results.minMidgetRGCToConeRatio, p.Results.maxMidgetRGCToConeRatio}, ...
        'StartPoint',p.Results.logitFitStartPoint, ...
        'Lower',[0,-2],'Upper',[10,0] );
    
    % Distribute the logitFit params to a loop variable
    fitParams(mm,1)=logitFit.slope;
    fitParams(mm,2)=logitFit.inflect;
    
    % Add the data to the figure
    if p.Results.makePlots
        plot(x,midgetRFtoConeRatio,p.Results.meridianSymbols{mm},'Color',[.7 .7 .7]);
        hold on
        xlabel('log10 proportion max cone density');
        ylabel('mRF density : cone density');
    end % if we are plotting
    
end % loop over meridians

% Calculate the median param values across meridians, but only for those
% meridians that we have declared hwe wish to use for this purpose
fitParams=median(fitParams(p.Results.meridiansIdxToUseForFitParams,:));

% Add a fit line to the plot
if p.Results.makePlots
    xFit= -2:.01:0;
    ylim([0 2.5]);
    plot( xFit, ...
        logisticFunc(fitParams(1), fitParams(2), p.Results.minMidgetRGCToConeRatio, p.Results.maxMidgetRGCToConeRatio, xFit),'-.m')
    legend({p.Results.meridianNames{:} 'fit'},'Location','southeast');
    title('midget RF : cone ratio as a function of relative cone density');
    drawnow
end


end % function


