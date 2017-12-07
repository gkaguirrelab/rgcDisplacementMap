function [ rgcDisplacementByMeridian, meridianAngleSupport, regularSupportPosDegRetina, opticDiscLocationByMeridian, mRF_RingCumulativeByMeridian, mRGC_RingCumulativeByMeridian, fitParamsByMeridian, fValsByMeridian, convergenceEccenDegRetinaByMeridian ] = createDisplacementModel( varargin )
% Creates a model of retinal ganglion cell displacement
%
% Description:
%   We begin with empirical measurements of cone and retinal ganglion cell
%   densities obtained from each of the four cardinal meridians of the
%   human retina. The data we use to build the model are from two papers
%   published by Curcio and colleagues in 1990.
%
%   We then engage in a modeling exercise to find a low-dimensional
%   parameterization that transforms cone density into midget ganglion
%   cell receptive field density (mRF) and retinal ganglion cell density
%   into midget retinal ganglion cell density (mRGC). In each case, the
%   models do not make use of explicit information regarding the retinal
%   position of the measurement to be transformed.
%
%   We are then in a position to model mRGC and mRF density as a function
%   of cone and RGC density, subject to a small number of parameters. Each
%   of these functions in turn can be expressed as a cumulative mRGC and
%   mRF count across the radial eccentricity of the retina.
%
%   As described by Drasdo (2007), the cumulative counts of mRF and mRGCs
%   will become equal at the eccentricity at which the RGCs are no longer
%   displaced. Further, the mRF and mRGC cumulative counts should be
%   equivalent beyond this point.
%
%   We engage in a non-linear model fit of the parameters that transform
%   cone density --> mRF density and RGC density --> mRGC density, subject
%   to the constraint that the cumulative mRF density must be greater than
%   the cumulative mRGC density within the displacement zone, and we
%   minimize an error function defined as the difference in mRGC and mRF
%   cumulative counts beyond the displacement zone. We set a convergence
%   target based upon the empirical measurements of Drasdo 2007.
%
%   We find that a single set of parameters that governs the RGC --> mRGC
%   transform is sufficient to model the data from all four meridians.
%   Further, we find sets of parameters that vary only slightly between the
%   meridians in the transform of cone density --> mRF density.
%
% Notes:
%   Units - All position values are on the retinal field, assuming a
%   spherical eye with the origin at the visual axis. Polar angle is
%   assigned zero for the nasal meridian on the retina, 90 degrees for the
%   superior meridian on the retina, etc. Eccentricity is defined in
%   degrees of distance from the visual axis along the spherical retina
%   away from the visual axis. The radius of the eye in mm is set as a
%   parameter in the unit conversion routines. Some routines internally
%   represent data in different units, including in degrees in visual field
%   and mm on the retina.
%
% Inputs:
%   None
%
% Optional key/value pairs:
%  'sampleResolutionDegreesRetina' - The calculations are performed across
%                           a regular sampling of eccentricity. This param
%                           defines sample resolution The sample resolution
%                           must be sufficient fine so that the cumulative
%                           is an accurate estimate of the integral.
%  'maxModeledEccentricityDegreesRetina' - The eccentricity extent of the
%                           model. This value must be sufficiently beyond
%                           the convergence zone so that there is a portion
%                           of the cumulative to match between the mRF and
%                           mRGC functions, but not so large as to venture
%                           into the periphery where our transform models
%                           are less accurate. We find that our results
%                           depend in unpredictable ways on the particular
%                           value selected, which is unfortunate.
%  'targetConvergenceOnCardinalMeridiansDegRetina' - The point in degrees
%                           at which displacement should become zero for
%                           each cadinal meridian. The default values were
%                           set by taking the apparent convergence points
%                           from Figure 2 of Drasdo 2007, and converting
%                           the values from mm retina to degrees retina.
%                           This provided values for the nasal and temporal
%                           meridians, and the inferior and superior
%                           meridian target values are simply the
%                           midpoints.
%  'targetMaxDisplacementDegRetina' - The desired maximum displacement
%                           value in retinal degrees. Taken from the maxium
%                           y-axis value from Figure 2 of Drasdo 2007,
%                           converted from mm retina to deg retina.
%  'meridianAngleResolutionDeg' - The resolution across polar angle for
%                           which displacements are calculated.
%  'displacementMapPixelsPerDegRetina' - The resolution in pixels per degree
%                           in which the displacement map is rendered.
%  'cone_to_mRF_linkTolerance' - We derive initial parameters for the
%                           linking function that converts cone density to
%                           mRF density. These values are used as the
%                           starting point for a search that aims to
%                           produce convergence of the mRF and mRGC
%                           cumulative functions (as well as other
%                           constraints). This parameter defines the upper
%                           and lower bounds for the search across the
%                           parameters as a multiplicative function.
%  'rgc_to_mRGC_linkTolerance' - Same as the prior param, except for the
%                           linking function of mRGC to RGC density.
%  'rgcLinkingFunctionFlavor' - Can use either the Dacey or Drasdo model to
%                           define the mapping of RGC --> mRGC
%  'rfInitialTransformParams' - Values for the linking function. If set to
%                           empty, these will be calculated using the
%                           develop model routine.
%  'rgcDrasdoInitialTransformParams' - Values for the linking function. If
%                           set to empty, these will be calculated using
%                           the develop model routine.
%  'rgcDaceyInitialTransformParams' - Values for the linking function. If
%                           set to empty, these will be calculated using
%                           the develop model routine. We find that our
%                           model requires parameters quite different from
%                           those that would be suggested by Dacey's
%                           values, so we define the default parameters
%                           here.
%  'coneDensityDataFileName' - The full path to a file (in standard form)
%                           the provides cone densities along the cardinal
%                           meridians. If left empty, the computed average
%                           Curcio values are used. Note that this setting
%                           controls the input to the calculation of
%                           displacement. The initial definition of the
%                           linking function parameters are not changed by
%                           this setting.
%  'rgcDensityDataFileName' - As above, but for RGC density.
%  'verbose'              - Do we give you the text?
%
%
% Outputs:
%   rgcDisplacementByMeridian - An m x p matrix, where m is number of 
%                           meridians and p is the number of eccentricity
%                           positions (starting from zero) modeled. The
%                           values are the displacement in retinal degrees
%                           of the RGC soma away from the fovea.
%   meridianAngleSupport  - An m x 1 vector of polar angle values (in
%                           degrees) that are the support for the meridian
%                           positions modeled.
%   regularSupportPosDegRetina - A 1 x p vector that contains the
%                           eccentricity in retinal degrees at which the
%                           model was evaluated along each meridian
%   mRF_cumulativeByMeridian - An m x p matrix. The cumulative RF counts
%                           across areal rings at each location.
%   mRGC_cumulativeByMeridian - An m x p matrix. The cumulative cell counts
%                           across areal rings at each location.
%   opticDiscLocationByMeridian - An m x p matrix. Is zero everywhere
%                           except at positions within the optic where it
%                           has a value of unity.
%   fitParamsByMeridian   - An m x k matrix, where m is the number of
%                           meridians and k is the number of parameters
%                           that define the linking functions that
%                           transform cone --> mRF and RGC --> mRGC.
%   fValsByMeridian       - An m x 1 vector with the values providing the
%                           error in minimizing the difference between mRF
%                           and mRGC cumulatives beyond the convergence
%                           point.
%   convergenceEccenDegRetinaByMeridian - An m x 1 vector with the values
%                           providing the location (in retinal degrees) at
%                           whicih the mRGC and mRF cumulative functions
%                           were found to converge.
%

%% Parse input and define variables
p = inputParser; p.KeepUnmatched = true;

% Optional anaysis params
p.addParameter('sampleResolutionDegreesRetina',0.01,@isnumeric);
p.addParameter('maxModeledEccentricityDegreesRetina',30,@isnumeric);
p.addParameter('targetConvergenceOnCardinalMeridiansDegRetina',[17 19.5 22 19.5],@isnumeric);
p.addParameter('targetMaxDisplacementDegRetina',3.2,@isnumeric);
p.addParameter('cardinalMeridianAngles',[0 90 180 270],@isnumeric);
p.addParameter('meridianAngleResolutionDeg',45,@isnumeric);
p.addParameter('displacementMapPixelsPerDegRetina',10,@isnumeric);
p.addParameter('cone_to_mRF_linkTolerance',1.1,@isnumeric);
p.addParameter('rgc_to_mRGC_linkTolerance',1.01,@isnumeric);
p.addParameter('rgcLinkingFunctionFlavor','Dacey',@(x)(stcmp(x,'Drasdo') | stcmp(x,'Dacey')));
p.addParameter('rfInitialTransformParams',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('rgcDrasdoInitialTransformParams',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('rgcDaceyInitialTransformParams',[4.2857 1.6],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('coneDensityDataFileName', [], @(x)(isempty(x) | ischar(x)));
p.addParameter('rgcDensityDataFileName', [], @(x)(isempty(x) | ischar(x)));

% Optional display params
p.addParameter('verbose',false,@islogical);

% parse
p.parse(varargin{:})


%% Setup
% Prepare the regular eccentricity support base
regularSupportPosDegRetina = ...
    0:p.Results.sampleResolutionDegreesRetina:p.Results.maxModeledEccentricityDegreesRetina;

% Prepare the set of meridian angles for which we will calculate
% displacement. This is a column vector.
meridianAngleSupport = (0:p.Results.meridianAngleResolutionDeg:(360-p.Results.meridianAngleResolutionDeg))';

% Prepare the vector targetConvergenceDegRetinaByMeridian. We interpolate
% from the values given for the cardinal merdians to the other meridians.
% We wrap the value at zero degrees around to 360 degrees to provide a
% circular interpolation.
% This is a column vector.
targetConvergenceByMeridianFit = fit(...
    [p.Results.cardinalMeridianAngles 360]',... % 
    [p.Results.targetConvergenceOnCardinalMeridiansDegRetina p.Results.targetConvergenceOnCardinalMeridiansDegRetina(1)]',...
    'pchipinterp');
targetConvergenceDegRetinaByMeridian = targetConvergenceByMeridianFit(meridianAngleSupport);

% Pre-allocated loop variables. We do not pre-allocate fitParamsByMeridian
% as the number of returned params could differ depending upon whether we
% go with Dacey or Drasdo model flavor
fValsByMeridian = zeros(length(meridianAngleSupport),1);
convergenceEccenDegRetinaByMeridian = zeros(length(meridianAngleSupport),1);
mRGC_RingCumulativeByMeridian = zeros(length(meridianAngleSupport),length(regularSupportPosDegRetina));
mRF_RingCumulativeByMeridian = zeros(length(meridianAngleSupport),length(regularSupportPosDegRetina));
rgcDisplacementByMeridian = zeros(length(meridianAngleSupport),length(regularSupportPosDegRetina));
opticDiscLocationByMeridian = zeros(length(meridianAngleSupport),length(regularSupportPosDegRetina));

%% Setup the initial values for the linking function parameters

% Derive parameters for the transformation of cone density to mRF density
if isempty(p.Results.rfInitialTransformParams)
    [ rfInitialTransformParams ] = developMidgetRFFractionModel();
else
    rfInitialTransformParams = p.Results.rfInitialTransformParams;
end

% Derive parameters for the transformation of RGC density to mRGC density
% Note that we can use either the Drasdo or Dacey base expression here
switch p.Results.rgcLinkingFunctionFlavor
    case 'Drasdo'
        % Derive linking function initial values if not passed
        if isempty(p.Results.rgcDrasdoInitialTransformParams)
            [ rgcInitialTransformParams ] = developDrasdoMidgetRGCFractionModel();
        else
            rgcInitialTransformParams = p.Results.rgcDrasdoInitialTransformParams;
        end
    case 'Dacey'
        % Derive linking function initial values if not passed
        if isempty(p.Results.rgcDaceyInitialTransformParams)
            [ rgcInitialTransformParams ] = developDaceyMidgetRGCFractionModel();
        else
            rgcInitialTransformParams = p.Results.rgcDaceyInitialTransformParams;
        end
    otherwise
        error('This is not an RGC linking function flavor that I know');
end


%% Loop over the meridians
for mm = 1:length(meridianAngleSupport)
    
    %% mRF_cumulative function
    % We build a function that returns the cumulative mRF density, subject
    % to two fit parameters. This function is based upon a model of cone
    % density.
    
    % Obtain a spline fit to the specified empirical density file
    [fitConeDensitySqDegRetina] = ...
        getSplineFitToConeDensitySqDegRetina(meridianAngleSupport(mm), ...
        'coneDensityDataFileName',p.Results.coneDensityDataFileName);
    % Define a variable with cone density over regular support and zero
    % values at the optic disc positions
    coneDensityOverRegularSupport = ...
        zeroOpticDiscPoints(fitConeDensitySqDegRetina(regularSupportPosDegRetina),regularSupportPosDegRetina, meridianAngleSupport(mm));
    % Create an anonymous function that returns mRF density as a function
    % of cone density, with the transform defined by the first two fitParams
    mRFDensityOverRegularSupport = ...
        @(fitParams) transformConeToMidgetRFDensity(coneDensityOverRegularSupport, ...
        'linkingFuncParams',fitParams(1:2));
    % Define anonymous function for the cumulative sum of mRF density
    mRF_RingCumulative = @(fitParams) calcRingCumulative(regularSupportPosDegRetina, mRFDensityOverRegularSupport(fitParams));
    
    
    %% mRGC_cumulative function
    % We build a function that returns the cumulative mRGC density, subject
    % to three fit parameters. This function is based upon a model of RGC
    % density.
    
    % Obtain a spline fit to the empirical RGC density data of Curcio 1990
    fitRGCDensitySqDegRetina = ...
        getSplineFitToRGCDensitySqDegRetina(meridianAngleSupport(mm), ...
        'rgcDensityDataFileName',p.Results.rgcDensityDataFileName);
    % Define a variable with RGC density over regular support and zero
    % values at the optic disc positions
    RGCDensityOverRegularSupport = ...
        zeroOpticDiscPoints(fitRGCDensitySqDegRetina(regularSupportPosDegRetina),regularSupportPosDegRetina, meridianAngleSupport(mm));

    % Create an anonymous function that returns mRGC density as a function
    % of RGC density, with the transform defined by the last two or three
    % fitParams This function can be set up with Drasdo or Dacey flavor
    switch p.Results.rgcLinkingFunctionFlavor
        case 'Drasdo'
            mRGCDensityOverRegularSupport = ...
                @(fitParams) transformRGCToMidgetRGCDensityDrasdo(regularSupportPosDegRetina,RGCDensityOverRegularSupport,...
                'linkingFuncParams',fitParams(3:end));
        case 'Dacey'
            mRGCDensityOverRegularSupport = ...
                @(fitParams) transformRGCToMidgetRGCDensityDacey(regularSupportPosDegRetina,RGCDensityOverRegularSupport,...
                'linkingFuncParams',fitParams(3:end));
        otherwise
            error('This is not an RGC linking function flavor that I know');
    end

    % Define anonymous function for the cumulative sum of mRGC density
    mRGC_RingCumulative = @(fitParams) calcRingCumulative(regularSupportPosDegRetina, mRGCDensityOverRegularSupport(fitParams));
    
    
    %% Non-linear constraint and error functions
    % Create a non-linear constraint that tests if the RF cumulative values
    % are greater than the RGC cumulative values at eccentricities less
    % than the displacement point
    nonlinconst = @(fitParams) constrainCumulativeAndDisplacement(regularSupportPosDegRetina, mRF_RingCumulative(fitParams), mRGC_RingCumulative(fitParams), targetConvergenceDegRetinaByMeridian(mm), p.Results.targetMaxDisplacementDegRetina);
    
    % The error function acts to minimize the diffrence between the
    % mRF and mRGC cumulative functions past the displacement point
    errorFunc = @(fitParams) errorMatchingRFandRGC(regularSupportPosDegRetina, mRF_RingCumulative(fitParams), mRGC_RingCumulative(fitParams), targetConvergenceDegRetinaByMeridian(mm));
    
    
    %% Perform the fit
    % We will search over the mRF and mRGC linking function parameters. We
    % set the upper and lower bounds as multipliers on the median values
    % found across meridians (with a bit of sign exponent trickery to
    % handle the direction of negative params)
    lb = [rfInitialTransformParams./(p.Results.cone_to_mRF_linkTolerance.^sign(rfInitialTransformParams)) rgcInitialTransformParams./(p.Results.rgc_to_mRGC_linkTolerance.^sign(rgcInitialTransformParams))];
    ub = [rfInitialTransformParams.*(p.Results.cone_to_mRF_linkTolerance.^sign(rfInitialTransformParams)) rgcInitialTransformParams.*(p.Results.rgc_to_mRGC_linkTolerance.^sign(rgcInitialTransformParams))];
    x0 = [rfInitialTransformParams rgcInitialTransformParams];
    
    % Set up the options
    options = optimoptions('fmincon', 'Algorithm','sqp','Display', 'none', 'ConstraintTolerance', 0.1);
    
    % Fit that sucker
    [fitParamsByMeridian(mm,:), fValsByMeridian(mm)]  = fmincon(errorFunc,x0,[],[],[],[],lb,ub,nonlinconst,options);
    
    % Calculate and store the cumulative, displacement, and optic disc fxns
    mRGC_RingCumulativeByMeridian(mm,:) = mRGC_RingCumulative(fitParamsByMeridian(mm,:));
    mRF_RingCumulativeByMeridian(mm,:) = mRF_RingCumulative(fitParamsByMeridian(mm,:));
    rgcDisplacementByMeridian(mm,:) = ...
        calcDisplacement(regularSupportPosDegRetina, mRF_RingCumulative(fitParamsByMeridian(mm,:)), mRGC_RingCumulative(fitParamsByMeridian(mm,:)));
    opticDiscLocationByMeridian(mm,findOpticDiscPositions(regularSupportPosDegRetina, meridianAngleSupport(mm)))=1;
    
    % Report the results for this meridian
    if p.Results.verbose
        zeroPoints = find(rgcDisplacementByMeridian(mm,:)==0);
        convergenceIdx = find(regularSupportPosDegRetina(zeroPoints) > 2,1);
        if isempty(convergenceIdx)
            outLine = ['Polar angle: ' num2str(meridianAngleSupport(mm)) ', max RGC displacement: ' num2str(max(rgcDisplacementByMeridian(mm,:))) ', target convergence: ' num2str(targetConvergenceDegRetinaByMeridian(mm)) ', found convergence: FAILED TO CONVERGE' ', error: ' num2str(round(fValsByMeridian(mm))) '\n'];
            fprintf(outLine);
        else
            convergenceEccenDegRetinaByMeridian(mm) = regularSupportPosDegRetina(zeroPoints(convergenceIdx));
            outLine = ['Polar angle: ' num2str(meridianAngleSupport(mm)) ', max RGC displacement: ' num2str(max(rgcDisplacementByMeridian(mm,:))) ', target convergence: ' num2str(targetConvergenceDegRetinaByMeridian(mm)) ', found convergence: ' num2str(convergenceEccenDegRetinaByMeridian(mm)) ', RMSE: ' num2str(round(fValsByMeridian(mm))) '\n'];
            fprintf(outLine);
        end
    end
    
end % loop over meridians


end % createDisplacementModel


%% LOCAL FUNCTIONS

function vectorOut = zeroOpticDiscPoints(vectorIn, regularSupportPosDegRetina, polarAngle)
	opticDiscIndices = findOpticDiscPositions(regularSupportPosDegRetina, polarAngle);
    vectorOut = vectorIn;
    vectorOut(opticDiscIndices) = 0;
end


function [c,ceq] = constrainCumulativeAndDisplacement(regularSupportPosDegRetina, mRF_RingCumulative, mRGC_RingCumulative, targetConvergencePointDegRetina, targetMaxDisplacementDegRetina)
% We implement two non-linear constraints. The first is that the mRGC
% cumulative values should not exceed the mRF cumulative values at any
% point prior to the convergence point. The second constraint is that the
% degree of displacement at any point should not exceed the maximum allowed
% displacement.

% First constraint:
% If there are any cumulative RGC values that are greater than the RF
% values at eccentricities less than the displacementPoint, then this
% violates the nonlinear fit constraint
withinRangeIdx = find(regularSupportPosDegRetina < targetConvergencePointDegRetina);
c = nansum(mRGC_RingCumulative(withinRangeIdx) > mRF_RingCumulative(withinRangeIdx));

% Second constraint:
% Calculate how many values in the displacement function exceed the maximum
% desired displacement.
displaceInDeg = calcDisplacement(regularSupportPosDegRetina, mRF_RingCumulative, mRGC_RingCumulative);
ceq = nansum(displaceInDeg > targetMaxDisplacementDegRetina);

end


function error = errorMatchingRFandRGC(regularSupportPosDeg, mRF_RingCumulative, mRGC_RingCumulative, targetConvergencePointDegRetina)
% The error is calculated as the RMSE of difference between the ring
% cumulative RF and RGC counts at retinal positions past the point where
% displacement should have ended

withinRangeIdx = find(regularSupportPosDeg > targetConvergencePointDegRetina);
error = sqrt(nanmean((mRGC_RingCumulative(withinRangeIdx) - mRF_RingCumulative(withinRangeIdx)).^2));
end


