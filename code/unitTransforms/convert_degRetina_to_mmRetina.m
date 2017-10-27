function supportPosMmRetina = convert_degRetina_to_mmRetina(supportPosDegRetina, varargin)
% convert_degRetina_to_mmRetina
%
%   Converts degrees on the retina to mm on the retina, assuming a
%   spherical eye. The value for the radius of curvature for the retina is
%   taken from the Table in Drasdo & Fowler 1974, for the value "Centre of
%   curvature of the retina"
%
% Input: 
%   supportPosDegRetina - Retinal postion(s) in degrees. Either scalar
%       value or vector input accepted.
%
% Output: 
%   supportPosMmRetina - Retinal postion(s) in degrees
%

%% Parse input and define variables
p = inputParser;

% Optional anaysis params
p.addParameter('radiusCurvatureRetinaMm',11.95,@isnumeric);

% parse
p.parse(varargin{:})

mmPerDegree = (2 * pi * p.Results.radiusCurvatureRetinaMm)/360;

supportPosMmRetina = supportPosDegRetina .* mmPerDegree;

end

