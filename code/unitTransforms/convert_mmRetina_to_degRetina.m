function supportPosDegRetina = convert_mmRetina_to_degRetina(supportPosMmRetina, varargin)
% Converts mm on the retina to degrees on the retina
%
% Description:
%   Unit conversion from mm on the retina to degrees across the retina.
%   The transformation assumes a spherical eye. The value used here is
%   taken from Drasdo & Fowler 1974.
%
% Inputs: 
%   supportPosMmRetina    - Retinal postion(s) in milimeters. Either scalar
%                           value or vector input accepted.
% Optional key/value pairs:
%  'radiusCurvatureRetinaMm'  - The value for the radius of curvature for 
%                           the retina. The Table in Drasdo & Fowler 1974
%                           offers a value of 11.95 for the "Centre of
%                           curvature of the retina". It appears that the
%                           Curcio 1990 papers use a value close to 11.5566
%
% Outputs: 
%   supportPosDegRetina   - Retinal postion(s) in degrees
%

%% Parse input and define variables
p = inputParser;

% Optional anaysis params
p.addParameter('radiusCurvatureRetinaMm',11.95,@isnumeric);

% parse
p.parse(varargin{:})

mmPerDegree = (2 * pi * p.Results.radiusCurvatureRetinaMm)/360;

supportPosDegRetina = supportPosMmRetina ./ mmPerDegree;

end

