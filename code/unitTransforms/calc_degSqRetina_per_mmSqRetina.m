function degSqRetinaPerMmSqRetina = calc_degSqRetina_per_mmSqRetina(varargin)
% Conversion factor for units of mm^2 to degrees^2 on the retina
%
% Description:
%   Returns the factor that converts degrees square on the retina to mm
%	square on the retina, assuming a spherical eye. The value for the
%	radius of curvature for the retina is taken from the Table in Drasdo &
%	Fowler 1974, for the value "Centre of curvature of the retina"
%
% Inputs:
%   none
%
% Optional key/value pairs:
%  'radiusCurvatureRetinaMm'  - The value for the radius of curvature for 
%                           the retina. The Table in Drasdo & Fowler 1974
%                           offers a value of 11.95 for the "Centre of
%                           curvature of the retina". It appears that the
%                           Curcio 1990 papers use a value close to 11.5566
%
% Outputs: 
%   degSqRetinaPerMmSqRetina - The conversion factor
%

%% Parse input and define variables
p = inputParser;

% Optional anaysis params
p.addParameter('radiusCurvatureRetinaMm',11.95,@isnumeric);

% parse
p.parse(varargin{:})

squareDegInASphere = 4 * pi * (360 / (2 * pi))^2;
squareMmInASphere = 4 * pi * p.Results.radiusCurvatureRetinaMm^2;

degSqRetinaPerMmSqRetina = squareDegInASphere / squareMmInASphere;

end

