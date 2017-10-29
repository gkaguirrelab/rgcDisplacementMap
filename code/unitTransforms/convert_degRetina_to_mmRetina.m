function supportPosMmRetina = convert_degRetina_to_mmRetina(supportPosDegRetina, varargin)
% convert_degRetina_to_mmRetina
%
%   Converts degrees on the retina to mm on the retina, assuming a
%   spherical eye.
%
% Input: 
%   supportPosDegRetina - Retinal postion(s) in degrees. Either scalar
%       value or vector input accepted.
%
% Output: 
%   supportPosMmRetina - Retinal postion(s) in degrees
%
% Optional:
%   radiusCurvatureRetinaMm - The value for the radius of curvature for the
%       retina. The Table in Drasdo & Fowler 1974 offers a value of 11.95
%       for the "Centre of curvature of the retina". It appears that the
%       Curcio 1990 papers use a value close to 11.5566

%% Parse input and define variables
p = inputParser;

% Optional anaysis params
p.addParameter('radiusCurvatureRetinaMm',11.5566,@isnumeric);


% parse
p.parse(varargin{:})

mmPerDegree = (2 * pi * p.Results.radiusCurvatureRetinaMm)/360;

supportPosMmRetina = supportPosDegRetina .* mmPerDegree;

end


