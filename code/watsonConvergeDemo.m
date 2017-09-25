function [ displacementMapDeg, meridianAngles, rgcDisplacementEachMeridian, mRGC_cumulativeEachMeridian, mRF_cumulativeEachMeridian ] = watsonConvergeDemo( varargin )
% watsonConvergeDemo( varargin )
%
% This is a modification of the primary routine (makeDisplacmenMap) that
% has as its purpose an illustration of the inability of the standard
% Watson equations to produce a well-structure displacement map
%
%
% OUTPUT
%   displacementMapDeg - The map of RGC radial displacement, in Cartesian
%       coordinates.
%   meridianAngles - a vector of polar angle values (in degrees) for which
%       the displacement values were calculated
%
% OPTIONS
%   sampleResolutionDegrees, maxModeledEccentricity - The calculations are
%       performed across a regular sampling of eccentricity. These params
%       deine sample resolution and the max modeled eccentricity. We
%       note that the sample resolution must be sufficient fine so that the
%       cumulative is an accurate estimate of the integral. Further, we
%       find that our results depend in unpredictable ways on the
%       particular maxModeledEccentricity selected. This latter value must
%       be sufficiently outside the displacement zone so that there is a
%       portion of the cumulative to match between the mRF and mRGC
%       functions, but not so large as to venture into the periphery where
%       our transform models areless accurate
%   targetDisplacementPointDeg - This is point in degrees at which
%       displacement should become zero for each cadinal meridian
%   meridianAngleResolutionDeg - The resolution across polar angle for
%       which displacements are calculated.
%   displacementMapPixelsPerDeg - The resolution in pixels per degree at
%       which the displacement map is rendered.
%   verbose - Do we give you the text?
%   makePlots - Do we make a figure?

%% Parse input and define variables
p = inputParser;

% Optional anaysis params
p.addParameter('sampleResolutionDegrees',0.01,@isnumeric);
p.addParameter('maxModeledEccentricity',30,@isnumeric);
p.addParameter('cardinalMeridianAngles',[0 90 180 270],@isnumeric);
p.addParameter('meridianAngleResolutionDeg',90,@isnumeric);
p.addParameter('displacementMapPixelsPerDeg',10,@isnumeric);

% Optional display params
p.addParameter('verbose',true,@islogical);

% parse
p.parse(varargin{:})


%% Setup
% Prepare the regular eccentricity support base
regularSupportPosDeg = ...
    0:p.Results.sampleResolutionDegrees:p.Results.maxModeledEccentricity;

% Prepare the set of meridian angles for which we will calculate
% displacement
meridianAngles = 0:p.Results.meridianAngleResolutionDeg:(360-p.Results.meridianAngleResolutionDeg);


%% Loop over the meridians
for mm = 1:length(meridianAngles)
    
    %% mRF_cumulative function
    % We build a function that returns the cumulative mRF density.
    
    % Start with Watson (2014) eq 8.
    mRFDensityOverRegularSupport = calcWatsonMidgetRFDensityByEccen(regularSupportPosDeg, meridianAngles(mm));
    % Define the cumulative sum of mRF density
    mRF_cumulative = calcCumulative(regularSupportPosDeg, mRFDensityOverRegularSupport);
    
    
    %% mRGC_cumulative function
    % We build a function that returns the cumulative mRGC density    
    % Obtain a spline fit to the empirical RGC density data of Curcio 1990
    RGCDensityFit = getSplineFitToRGCDensity(meridianAngles(mm));
    % Create a function that returns mRGC density as a function of
    % RGC density, with the transform defined by Watson (2014) Eq 7 (which
    % itself is taken from Drasdo 2017)
    mRGCDensityOverRegularSupport = ...
        RGCDensityFit(regularSupportPosDeg)' .* ...
        calcWatsonMidgetFractionByEccen(regularSupportPosDeg,0.8928,41.03);
    % Define the cumulative sum of mRGC density
    mRGC_cumulative = calcCumulative(regularSupportPosDeg, mRGCDensityOverRegularSupport);
        
    % Calculate and store displacement and cumulative functions
    mRGC_cumulativeEachMeridian(mm,:)=mRGC_cumulative;
    mRF_cumulativeEachMeridian(mm,:)=mRF_cumulative;
    rgcDisplacementEachMeridian(mm,:)=calcDisplacement(regularSupportPosDeg, mRGC_cumulative, mRF_cumulative);
    
    % Report the results for this meridian
    if p.Results.verbose
        zeroPoints = find(rgcDisplacementEachMeridian(mm,:)==0);
        convergenceIdx = find(regularSupportPosDeg(zeroPoints) > 2,1);
        convergenceEccen = regularSupportPosDeg(zeroPoints(convergenceIdx));
        outLine = ['Polar angle: ' num2str(meridianAngles(mm)) ', max RGC displacement: ' num2str(max(rgcDisplacementEachMeridian(mm,:))) ', found convergence: ' num2str(convergenceEccen) '\n'];
        fprintf(outLine);
    end
    
end % loop over meridians

% Create the displacement map
imRdim = p.Results.maxModeledEccentricity * p.Results.displacementMapPixelsPerDeg * 2;
maxDisplacementDeg = max(rgcDisplacementEachMeridian(:));
imP=rgcDisplacementEachMeridian'./maxDisplacementDeg;
imR = PolarToIm (imP, 0, 1, imRdim, imRdim);
displacementMapDeg = imrotate(imR .* maxDisplacementDeg,-90);



end % watsonConvergeDemo


%% LOCAL FUNCTIONS

function displaceInDeg = calcDisplacement(regularSupportPosDeg, countPerRingRGC, countPerRingRF)

% Determine the sample resolution by a difference operation
tmp = diff(regularSupportPosDeg);
sampleResolutionDegrees = tmp(1);
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
displaceInDeg = ((1:1:length(displaceInSamples))-displaceInSamples ) * sampleResolutionDegrees;
% Zero out negative values after the point of convergence
displaceInDeg(find(displaceInDeg < 0,1):end)=0;

end

