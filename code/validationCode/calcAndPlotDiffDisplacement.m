function [rgcDiffDisplacementEachMeridian,diffDisplacementMapDeg] = calcAndPlotDiffDisplacement(rgcDisplacementEachMeridian,varargin)



%% Parse input and define variables
p = inputParser;

% Optional anaysis params
p.addParameter('maxModeledEccentricity',30,@isnumeric);
p.addParameter('displacementMapPixelsPerDeg',10,@isnumeric);
p.addParameter('verbose',true,@islogical);

% parse
p.parse(varargin{:})


%% Check neighbor diff

rgcDiffDisplacementEachMeridian =  diff(rgcDisplacementEachMeridian,1,2);

if p.Results.verbose
    figure
    imagesc(rgcDiffDisplacementEachMeridian);
end

% Convert to radial image
imRdim = p.Results.maxModeledEccentricity *  p.Results.displacementMapPixelsPerDeg * 2;
maxDiffDisplacementDeg = max(rgcDiffDisplacementEachMeridian(:));
imP=rgcDiffDisplacementEachMeridian'./maxDiffDisplacementDeg;
imR = PolarToIm (imP, 0, 1, imRdim, imRdim);
diffDisplacementMapDeg = imrotate(imR .* maxDiffDisplacementDeg,-90);

if p.Results.verbose
    figure
    imagesc(diffDisplacementMapDeg)
end