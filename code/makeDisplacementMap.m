function [ displacementMapDeg, fitParams, meridianAngles, rgcDisplacementEachMeridian, mRGC_cumulativeEachMeridian, mRF_cumulativeEachMeridian, convergenceEccen] = makeDisplacementMap( varargin )
% makeDisplacementMap -  This routine models retinal ganglion cell displacement.
%
% Our strategy is to begin with empirical measurements of cone and retinal
% ganglion cell densities obtained from each of the four cardinal meridians
% of the human retina. The data we use are from two papers published by
% Curcio and colleagues in 1990.
%
% We then engage in a modeling exercise to find a low-dimensional
% parameterization of the transformation of cone density into midget
% ganglion cell receptive field density (mRF) and of retinal ganglion cell
% density into midget retinal ganglion cell density (mRGC). In each case,
% the models do not make use of explicit information regarding the retinal
% position of the measurement to be transformed.
%
% We are then in a position to model mRGC and mRF density as a function of
% cone and RGC density, subject to a small number of parameters. Each of
% these functions in turn can be expressed as a cumulative mRGC and mRF
% count across the radial eccentricity of the retina.
%
% As described by Drasdo (2007), the cumulative counts of mRF and mRGCs
% will become equal at the eccentricity at which the RGCs are no longer
% displaced. Further, the mRF and mRGC cumulative counts should be
% equivalent beyond this point.
%
% We engage in a non-linear model fit of the parameters that transform
% cone density --> mRF density and RGC density --> mRGC density, subject
% to the constraint that the cumulative mRF density must be greater than
% the cumulative mRGC density within the displacement zone, and we minimize
% an error function defined as the difference in mRGC and mRF cumulative
% counts beyond the displacement zone. We set the displacement zone target
% to be 17 degrees eccentricity for all meridians, except for the nasal
% meridian within which the displacement point of 11 degrees is enforced by
% the presence of the optic nerve head.
%
% We find that a single set of parameters that governs the RGC --> mRGC
% transform is sufficient to model the data from all four meridians.
% Further, we find sets of parameters that vary only slightly between the
% meridians in the transform of cone density --> mRF density.
%
%
% OUTPUT
%   displacementMapDeg - The map of RGC radial displacement, in Cartesian
%       coordinates.
%   fitParams - The values of the five parameters that adjust the transform
%      of cone --> mRF and RGC --> mRGC, provided for each meridian
%   meridianAngles - a vector of polar angle values (in degrees) for which
%       the displacement values were calculated
%
% OPTIONS
%   sampleResolutionDegreesRetina - The calculations are
%       performed across a regular sampling of eccentricity. This param
%       defines sample resolution The sample resolution must be sufficient
%       fine so that the cumulative is an accurate estimate of the
%       integral.
%   maxModeledEccentricityRetina - The eccentricity extent of the model. This
%       value must be sufficiently outside the displacement zone so that
%       there is a portion of the cumulative to match between the mRF and
%       mRGC functions, but not so large as to venture into the periphery
%       where our transform models areless accurate. We find that our
%       results depend in unpredictable ways on the particular
%       maxModeledEccentricityRetina selected, which is unfortunate.
%   targetConvergenceOnCardinalMeridiansDegRetina - This is point in
%       degrees at which displacement should become zero for each cadinal
%       meridian. The default values were set by taking the apparent
%       convergence points from Figure 2 of Drasdo 2007, and converting
%       the values from mm retina to degrees retina. This provided values
%       for the nasal and temporal meridians, and the inferior and superior
%       meridian target values are simply the midpoints.
%   targetMaxDisplacementDegRetina - The desired maximum displacement
%       value in retinal degrees. Taken from the maxium y-axis value from
%       Figure 2 of Drasdo 2007, converted from mm retina to deg retina.
%   meridianAngleResolutionDeg - The resolution across polar angle for
%       which displacements are calculated.
%   displacementMapPixelsPerDegRetina - The resolution in pixels per degree at
%       which the displacement map is rendered.
%   cone_to_mRF_linkTolerance - We derive initial parameters for the
%       linking function that converts cone density to mRF density. These
%       values are used as the starting point for a search that aims to
%       produce convergence of the mRF and mRGC cumulative functions (as
%       well as other constraints). This parameter defines the upper and
%       lower bounds for the search across the parameters as a
%       multiplicative function.
%   rgc_to_mRGC_linkTolerance - Same as the prior param, except for the
%       linking function of mRGC to RGC density.
%   rfInitialTransformParams - values for the linking function. If set to
%       empty, these will be calculated using the develop model routine.
%   rgcDrasdoInitialTransformParams - values for the linking function. If set to
%       empty, these will be calculated using the develop model routine.
%
%   verbose - Do we give you the text?

%% Parse input and define variables
p = inputParser; p.KeepUnmatched = true;

% Optional anaysis params
p.addParameter('sampleResolutionDegreesRetina',0.01,@isnumeric);
p.addParameter('maxModeledEccentricityRetina',30,@isnumeric);
p.addParameter('targetConvergenceOnCardinalMeridiansDegRetina',[17 19.5 22 19.5],@isnumeric);
p.addParameter('targetMaxDisplacementDegRetina',3.2,@isnumeric);

p.addParameter('cardinalMeridianAngles',[0 90 180 270],@isnumeric);
p.addParameter('meridianAngleResolutionDeg',15,@isnumeric);
p.addParameter('displacementMapPixelsPerDegRetina',10,@isnumeric);
p.addParameter('cone_to_mRF_linkTolerance',1.2,@isnumeric);
p.addParameter('rgc_to_mRGC_linkTolerance',1.2,@isnumeric);
p.addParameter('rgcLinkingFunctionFlavor','Dacey',@(x)(stcmp(x,'Drasdo') | stcmp(x,'Dacey')));
p.addParameter('rfInitialTransformParams',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('rgcDrasdoInitialTransformParams',[],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('rgcDaceyInitialTransformParams',[5 1.35],@(x)(isempty(x) | isnumeric(x)));


% Optional display params
p.addParameter('verbose',false,@islogical);

% parse
p.parse(varargin{:})


%% Setup
% Prepare the regular eccentricity support base
regularSupportPosDegRetina = ...
    0:p.Results.sampleResolutionDegreesRetina:p.Results.maxModeledEccentricityRetina;

% Prepare the set of meridian angles for which we will calculate
% displacement
meridianAngles = 0:p.Results.meridianAngleResolutionDeg:(360-p.Results.meridianAngleResolutionDeg);

% Prepare the vector of targetDisplacement values. We interpolate from the
% values given for the cardinal merdians to the other meridians
targetValues = [p.Results.targetConvergenceOnCardinalMeridiansDegRetina p.Results.targetConvergenceOnCardinalMeridiansDegRetina(1)];
angleBase = [p.Results.cardinalMeridianAngles 360];
targetByAngleFit = fit(angleBase',targetValues','pchipinterp');
targetConvergenceDegRetinaByMeridian = targetByAngleFit(meridianAngles);


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
for mm = 1:length(meridianAngles)
    
    %% mRF_cumulative function
    % We build a function that returns the cumulative mRF density, subject
    % to two fit parameters. This function is based upon a model of cone
    % density.
    
    % Obtain a spline fit to the empirical cone density data of Curcio 1990
    [coneDensityFit] = getSplineFitToConeDensity(meridianAngles(mm));
    % Create an anonymous function that returns mRF density as a function
    % of cone density, with the transform defined by the first two fitParams
    mRFDensityOverRegularSupport = ...
        @(fitParams) transformConeToMidgetRFDensity(coneDensityFit(regularSupportPosDegRetina), ...
        'linkingFuncParams',fitParams(1:2))';
    % Define anonymous function for the cumulative sum of mRF density
    mRF_cumulative = @(fitParams) calcCumulative(regularSupportPosDegRetina, mRFDensityOverRegularSupport(fitParams));
    
    
    %% mRGC_cumulative function
    % We build a function that returns the cumulative mRGC density, subject
    % to three fit parameters. This function is based upon a model of RGC
    % density.
    
    % Obtain a spline fit to the empirical RGC density data of Curcio 1990
    rgcDensityFit = getSplineFitToRGCDensity(meridianAngles(mm));

    % Create an anonymous function that returns mRGC density as a function of
    % RGC density, with the transform defined by the last two or three fitParams
    % This function can be set up with Drasdo or Dacey flavor
    switch p.Results.rgcLinkingFunctionFlavor
        case 'Drasdo'
            mRGCDensityOverRegularSupport = ...
                @(fitParams) transformRGCToMidgetRGCDensityDrasdo(regularSupportPosDegRetina,rgcDensityFit(regularSupportPosDegRetina)',...
                'linkingFuncParams',fitParams(3:end));
        case 'Dacey'
            mRGCDensityOverRegularSupport = ...
                @(fitParams) transformRGCToMidgetRGCDensityDacey(regularSupportPosDegRetina,rgcDensityFit(regularSupportPosDegRetina)',...
                'linkingFuncParams',fitParams(3:end));
        otherwise
            error('This is not an RGC linking function flavor that I know');
    end

    % Define anonymous function for the cumulative sum of mRGC density
    mRGC_cumulative = @(fitParams) calcCumulative(regularSupportPosDegRetina, mRGCDensityOverRegularSupport(fitParams));
    
    
    %% Non-linear constraint and error functions
    % Create a non-linear constraint that tests if the RF cumulative values
    % are greater than the RGC cumulative values at eccentricities less
    % than the displacement point
    nonlinconst = @(fitParams) constrainCumulativeAndDisplacement(regularSupportPosDegRetina, mRF_cumulative(fitParams), mRGC_cumulative(fitParams), targetConvergenceDegRetinaByMeridian(mm), p.Results.targetMaxDisplacementDegRetina);
    
    % The error function acts to minimize the diffrence between the
    % mRF and mRGC cumulative functions past the displacement point
    errorFunc = @(fitParams) errorMatchingRFandRGC(regularSupportPosDegRetina, mRF_cumulative(fitParams), mRGC_cumulative(fitParams), targetConvergenceDegRetinaByMeridian(mm));
    
    
    %% Perform the fit
    % We will search over the mRF and mRGC linking function parameters. We
    % set the upper and lower bounds as multipliers on the median values
    % found across meridians (with a bit of sign exponent trickery to
    % handle the direction of negative params)
    lb = [rfInitialTransformParams./(p.Results.cone_to_mRF_linkTolerance.^sign(rfInitialTransformParams)) rgcInitialTransformParams./(p.Results.rgc_to_mRGC_linkTolerance.^sign(rgcInitialTransformParams))];
    ub = [rfInitialTransformParams.*(p.Results.cone_to_mRF_linkTolerance.^sign(rfInitialTransformParams)) rgcInitialTransformParams.*(p.Results.rgc_to_mRGC_linkTolerance.^sign(rgcInitialTransformParams))];
    x0 = [rfInitialTransformParams rgcInitialTransformParams];
    
    % Set up the options
    options = optimoptions('fmincon', 'Display', 'none', 'ConstraintTolerance', 0.1);
    
    % Fit that sucker
    fitParams(mm,:) = fmincon(errorFunc,x0,[],[],[],[],lb,ub,nonlinconst,options);
    
    % Calculate and store displacement and cumulative functions
    mRGC_cumulativeEachMeridian(mm,:)=mRGC_cumulative(fitParams(mm,:));
    mRF_cumulativeEachMeridian(mm,:)=mRF_cumulative(fitParams(mm,:));
    rgcDisplacementEachMeridian(mm,:)=calcDisplacement(regularSupportPosDegRetina, mRGC_cumulative(fitParams(mm,:)), mRF_cumulative(fitParams(mm,:)));
    
    % Report the results for this meridian
    if p.Results.verbose
        zeroPoints = find(rgcDisplacementEachMeridian(mm,:)==0);
        convergenceIdx = find(regularSupportPosDegRetina(zeroPoints) > 2,1);
        if isempty(convergenceIdx)
            outLine = ['Polar angle: ' num2str(meridianAngles(mm)) ', max RGC displacement: ' num2str(max(rgcDisplacementEachMeridian(mm,:))) ', target convergence: ' num2str(targetConvergenceDegRetinaByMeridian(mm)) ', found convergence: FAILED TO CONVERGE\n'];
            fprintf(outLine);
        else
            convergenceEccen(mm) = regularSupportPosDegRetina(zeroPoints(convergenceIdx));
            outLine = ['Polar angle: ' num2str(meridianAngles(mm)) ', max RGC displacement: ' num2str(max(rgcDisplacementEachMeridian(mm,:))) ', target convergence: ' num2str(targetConvergenceDegRetinaByMeridian(mm)) ', found convergence: ' num2str(convergenceEccen(mm)) '\n'];
            fprintf(outLine);
        end
    end
    
end % loop over meridians

% Create the displacement map
imRdim = (p.Results.maxModeledEccentricityRetina * p.Results.displacementMapPixelsPerDegRetina * 2)-1;
displacementMapDeg = convertPolarMapToImageMap(rgcDisplacementEachMeridian, imRdim);


end % calcDisplacementMap


%% LOCAL FUNCTIONS

function [c,ceq] = constrainCumulativeAndDisplacement(regularSupportPosDeg, countPerRingRF, countPerRingRGC, targetConvergencePointDegRetina, targetMaxDisplacementDegRetina)

% If there are any cumulative RGC values that are greater than the RF
% values at eccentricities less than the displacementPoint, then this
% violates the nonlinear fit constraint
withinRangeIdx = find(regularSupportPosDeg < targetConvergencePointDegRetina);
c = sum(countPerRingRGC(withinRangeIdx) > countPerRingRF(withinRangeIdx));

% Calculate how many values in the displacement function exceed the maximum
% desired displacement.
displaceInDeg = calcDisplacement(regularSupportPosDeg, countPerRingRGC, countPerRingRF);
ceq = sum(displaceInDeg > targetMaxDisplacementDegRetina);

end

function error = errorMatchingRFandRGC(regularSupportPosDeg, countPerRingRF, countPerRingRGC, targetConvergencePointDegRetina)

% The error is calculated as the SSQ of the absolute difference between the
% cumulative RF and RGC counts at retinal positions past the point where
% displacement should have ended

withinRangeIdx = find(regularSupportPosDeg > targetConvergencePointDegRetina);
error = sqrt(sum((countPerRingRGC(withinRangeIdx) - countPerRingRF(withinRangeIdx)).^2));
end

function displaceInDeg = calcDisplacement(regularSupportPosDeg, countPerRingRGC, countPerRingRF)

% Determine the sample resolution by a difference operation
tmp = diff(regularSupportPosDeg);
sampleResolutionDegreesRetina = tmp(1);
% Measure the displacement (in degrees). First, for each cumulative RGC
% density value, identify the array index of the first
% cumulative RF density value that is equal or greater.
displaceInSamples=arrayfun(@(x) find(countPerRingRF>=x,1), countPerRingRGC,'UniformOutput',false);
% Now some array operations to get these values out of cells and in to a
% numeric vector
emptyCells = find(cellfun(@(x) isempty(x), displaceInSamples));
displaceInSamples(emptyCells)={NaN};
displaceInSamples=cell2mat(displaceInSamples(cellfun(@(x) ~isempty(x), displaceInSamples)));
% The displacement in array samples is the array index of each cumulative
% RGC density measurement, minus the array index of the matching RF density
% measurement. This difference is then multiplied by the sample resolution
% to get the displacement in degrees.
displaceInDeg = ((1:1:length(displaceInSamples))-displaceInSamples ) * sampleResolutionDegreesRetina;
% Zero out negative values after the point of convergence
displaceInDeg(find(displaceInDeg < 0,1):end)=0;

end

