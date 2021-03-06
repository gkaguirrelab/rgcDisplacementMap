function [ fitParams, figHandle ] = developDaceyMidgetRGCFractionModel( varargin )
% Find parameters that transform RGC density to midget RGC density (Dacey)
%
% Syntax:
%  [ fitParams, figHandle ] = developDaceyMidgetRGCFractionModel( varargin )
%
% Description:
%   Our goal is to be able to derive the fraction of RGCs that are midget
%   RGCs at a point in the retina, but without reference to the polar angle
%   or eccentricity of the measurement.
%
%   To solve this problem, we first take the cumulative of the RGC density
%   function along a meridian. We then express this function as a
%   proportion relative to the total number of RGCs in the retina (1e6). 
%
%   Next, the proportional, cumulative RGC density function is related to
%   the predicted midget fraction for that location, as derived by a fit to
%   the data of Dacey 1993, Figure 19B. We observe that logisitic function
%   fits this relationship well. We lock the parameters for the maximum and
%   minimum midget fractions, leaving two free parameters for the fit.
%
%   We obtain the fit parameters for the midget fraction expression for
%   each of the four meridians for which we have data, and take the median
%   fit parameters across merdians. This gives us a functional form that
%   can be used to convert RGC density to midget RGC density knowing only
%   the proportion of the cumulative RGC density function at the location
%   for which the measurement was made. If we term this proportion value x,
%   then the relationship is:
%
%       midgetFraction = minMidgetFractionRatio + 
%           (maxMidgetFractionRatio - minMidgetFractionRatio) ./ 
%           (1+sign(x./inflect).*abs((x./inflect).^slope))
%
%   where inflect and slope are free paramters of the logisitic fit, and
%   maxMidgetFractionRatio and minMidgetFractionRatio are the values of the
%   midget fraction at the fovea and periphery, respectivey.
%
%   We have adopted a value of maxMidgetFractionRatio based upon the
%   measurements reported by Don Miller's group (Liu et al 2017 PNAS) that
%   show a peak midget fraction of ~0.95 close to the fovea. The minimum
%   midget fraction is taken from the asymptotic value of a fit to Dacey's
%   measurements.
%
% Inputs:
%   none
%
% Optional key/value pairs:
%  'referenceEccenDegVisual' - The reference eccentricity for the proportion 
%                           of the cumulative RGC density. The proportion
%                           function will have a value of unity at this
%                           point. We use 10� here for the practical reason
%                           that this is within the range for which we
%                           have OCT measurements of the RGC layer
%                           thickness, and we wish in the future to model
%                           such data using these functions.
%  'supportResolutionDegVisual' - the resolution (in degrees) for which
%                           the calculations are made.
%  'supportEccenMaxDegVisual' - The maximum eccentricity used for
%                           modeling
%  'meridianNames'        - Cell array of the text string names of the
%                           meridians
%  'meridianAngles'       - Polar angle values assigned to the meridians
%  'meridiansIdxToUseForFitParams' - Index of the meridians from which we 
%                           will calculate the median fit param values.
%  'maxMidgetFractionRatio' - The midget fraction assigned to the fovea
%  'minMidgetFractionRatio' - The midget fraction assigned to the far
%                           periphery
%  'logitFitStartPoint'   - Initial values used for the fit of the
%                           logisitic function.
%  'makePlots'            - Do we make a figure?
%
% Outputs:
%   fitParams             - The median values of [slope inflect] across the 
%                           merdians selected to be used
%   figHandle             - Handle to a figure showing the fits. Empty if 
%                           no plotting requested.
%
% Examples:
%{
    fitParams = developDaceyMidgetRGCFractionModel('makePlots',true)
%}
%{
    figure
    regularSupportPosDegVisual = 0:0.1:30;
    [~, daceyMidgetFraction, daceyDataSupportPosDegVisual] = calcDaceyMidgetFractionByEccenDegVisual(regularSupportPosDegVisual);
    plot(daceyDataSupportPosDegVisual,daceyMidgetFraction,'or');
    hold on
    fitParams = developDaceyMidgetRGCFractionModel();
    fitRGCDensitySqDegVisual = getSplineFitToRGCDensitySqDegVisual(180);
    [ ~, midgetFraction ] = transformRGCToMidgetRGCDensityDacey( regularSupportPosDegVisual, fitRGCDensitySqDegVisual(regularSupportPosDegVisual), 'linkingFuncParams', fitParams );
    plot(regularSupportPosDegVisual,midgetFraction)
%}

%% Parse input and define variables
p = inputParser;

% Optional anaysis params
p.addParameter('totalRetinalRGCs',1.5e6,@isscalar);
p.addParameter('supportResolutionDegVisual',0.01,@isnumeric);
p.addParameter('supportEccenMaxDegVisual',60,@isnumeric);
p.addParameter('meridianNames',{'Nasal' 'Superior' 'Temporal' 'Inferior'},@iscell);
p.addParameter('meridianAngles',[0, 90, 180, 270],@isnumeric);
p.addParameter('meridianSymbols',{'.','x','o','^'},@cell);
p.addParameter('meridiansIdxToUseForFitParams',[1 2 3 4],@isnumeric);
p.addParameter('minMidgetFractionRatio',0.41,@isnumeric);
p.addParameter('maxMidgetFractionRatio',0.95,@isnumeric);
p.addParameter('logitFitStartPoint',[5 0.5],@isnumeric);

% Optional display params
p.addParameter('makePlots',false,@islogical);

% parse
p.parse(varargin{:})


%% House keeping and setup

% Define a four-parameter logistic function that will be used to fit the
% modeled relationship. Two of the parameters (max and min asymptote) are
% locked by the passed parameter
logisticFunc = fittype( @(slope,inflect,minMidgetFractionRatio,maxMidgetFractionRatio,x) minMidgetFractionRatio+(maxMidgetFractionRatio-minMidgetFractionRatio)./(1+sign(x./inflect).*abs((x./inflect).^slope)), ...
    'independent','x','dependent','y','problem',{'minMidgetFractionRatio','maxMidgetFractionRatio'});

% Define the regular-spaced eccentricity support over which we will model
% the anatomical retinal functions
regularSupportPosDegVisual = ...
    0:p.Results.supportResolutionDegVisual:p.Results.supportEccenMaxDegVisual;

% Prepare a figure if requested
if p.Results.makePlots
    figHandle=figure;
else
    figHandle=[];
end

%% Loop over the meridians

for mm = 1:length(p.Results.meridianAngles)
    
	% Get a spline to the Curcio RGC density data for this meridian 
    fitRGCDensitySqDegVisual = getSplineFitToRGCDensitySqDegVisual(p.Results.meridianAngles(mm));

    % Obtain the ring cumulative RGC function
    RGC_ringcount = calcRingCumulative(regularSupportPosDegVisual,fitRGCDensitySqDegVisual(regularSupportPosDegVisual));
    
    % Calculate a proportion of the cumulative RGC density counts, relative
    % to the reference point (which is assigned a value of unity)
    propRGC_ringcount=RGC_ringcount./p.Results.totalRetinalRGCs;
    
    % Obtain the Dacey midget fraction as a function of eccentricity
    midgetFractionByEccenDegVisual = calcDaceyMidgetFractionByEccenDegVisual(regularSupportPosDegVisual)';
    
    % Adjust the Dacey fraction so that it has the specified max and
    % minimum values
    midgetFractionByEccenDegVisual = midgetFractionByEccenDegVisual-min(midgetFractionByEccenDegVisual);
    midgetFractionByEccenDegVisual = midgetFractionByEccenDegVisual./max(midgetFractionByEccenDegVisual);
    midgetFractionByEccenDegVisual = midgetFractionByEccenDegVisual.*(p.Results.maxMidgetFractionRatio-p.Results.minMidgetFractionRatio);
    midgetFractionByEccenDegVisual = midgetFractionByEccenDegVisual+p.Results.minMidgetFractionRatio;
    
    % Perform the fit and save the param values.
    % The x values proportionRGC
    % The y values are the midget fraction
    
    % turn off some warnings that are produced by the fit and plot
    warningState = warning;
    warning('off','curvefit:fit:complexXusingOnlyReal');
    warning('off','MATLAB:plot:IgnoreImaginaryXYPart');
    
    % Perform the logistic fit. Note that the max and min asymptote are
    % pinned by the passed parameters
    logitFit = fit(propRGC_ringcount',midgetFractionByEccenDegVisual',logisticFunc, ...
        'problem',{p.Results.minMidgetFractionRatio, p.Results.maxMidgetFractionRatio}, ...
        'StartPoint',p.Results.logitFitStartPoint, ...
        'Lower',[0,0],'Upper',[50,5] );
    
    % Distribute the logitFit params to a loop variable
    fitParams(mm,1)=logitFit.slope;
    fitParams(mm,2)=logitFit.inflect;
    
    % Add the data to the figure
    if p.Results.makePlots
        plot(propRGC_ringcount,midgetFractionByEccenDegVisual,p.Results.meridianSymbols{mm},'color',[.8 .8 .8])
        hold on
        ylim([0.4 1]);
        xlim([0 2]);
        xlabel('proportion cumulative RGC density count at reference eccen');
        ylabel('midget fraction');
    end % if we are plotting
    
    % restore the saved warning state
    warning(warningState);
    
end % loop over meridians

% Calculate the median param values across meridians, but only for those
% meridians that we have declared hwe wish to use for this purpose
fitParams=median(fitParams(p.Results.meridiansIdxToUseForFitParams,:));

% Add a fit line to the plot
if p.Results.makePlots
    xFit=0:.01:1;
    plot( xFit, ...
        logisticFunc(fitParams(1), fitParams(2), p.Results.minMidgetFractionRatio, p.Results.maxMidgetFractionRatio, xFit),'-.m')
    legend({p.Results.meridianNames{:} 'fit'},'Location','southwest');
    title('midget fraction as a function of relative RGC cumulative density');
    drawnow
end


end % function


