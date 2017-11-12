function degSqRetinaPerMmSqRetina = calc_degSqRetina_per_mmSqRetina(varargin)
% convert_mmRetina_to_degRetina
%
%   Returns the degrees square on the retina per mm square on the retina,
%   assuming a spherical eye. The value for the radius of curvature for the
%   retina is taken from the Table in Drasdo & Fowler 1974, for the value
%   "Centre of curvature of the retina"
%
% Output: 
%   degSqRetinaPerMmSqRetina - The conversion factor
%
% Optional
%
%   radiusCurvatureRetinaMm - The value for the radius of curvature for the
%       retina. The Table in Drasdo & Fowler 1974 offers a value of 11.95
%       for the "Centre of curvature of the retina". It appears that the
%       Curcio 1990 papers use a value close to 11.5566
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

