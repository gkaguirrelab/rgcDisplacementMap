function supportPosDegRetina = convert_mmRetina_to_degRetina(supportPosMmRetina, varargin)
% convert_mmRetina_to_degRetina
%
%   Converts mm on the retina to degrees on the retina, assuming a
%   spherical eye. The value for the radius of curvature for the retina is
%   taken from the Table in Drasdo & Fowler 1974, for the value "Centre of
%   curvature of the retina"
%
% Input: 
%   supportPosMmRetina  = Retinal postion(s) in milimeters. Either scalar value or
%                     vector input accepted.
%
% Output: 
%   supportPosDegRetina = Retinal postion(s) in degrees
%

%% Parse input and define variables
p = inputParser;

% Optional anaysis params
p.addParameter('radiusCurvatureRetinaMm',11.95,@isnumeric);

% parse
p.parse(varargin{:})

mmPerDegree = (p.Results.radiusCurvatureRetinaMm * 2 * pi)/360;

supportPosDegRetina = supportPosMmRetina .* mmPerDegree;

end

