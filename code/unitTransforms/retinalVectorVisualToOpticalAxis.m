function [ angleVisualToOpticalAxis, distanceMmRetinaVisualToOpticalAxis ] = retinalVectorVisualToOpticalAxis( varargin )
% retinalVectorVisualToOpticalAxis
%
% Returns the polar angle and distance in mm of the vector that arises at
% the visual axis of the eye and terminates at the optical axis of the eye.
%
% The coordinate system for polar angle is in degrees, and corresponds to
%   0=nasal;90=superior;180=temporal;270=inferior
% on the retinal field.
%
% Optional arguments:
%   displacementOpticalAxisFromVisualAlongNasalMeridianMm
%   displacementOpticalAxisFromVisualAlongSuperiorMeridianMm - From Watson
%       2014 (in turn from Charman 1991, quoting Emsley 1952), the optical
%       axis intersects the retina 1.5 mm nasal and 0.5 mm superior to the
%       visual axis.
%


%% Parse input and define variables
p = inputParser;

% Optional analysis params
p.addParameter('displacementOpticalAxisFromVisualAlongNasalMeridianMm',1.5,@isnumeric);
p.addParameter('displacementOpticalAxisFromVisualAlongSuperiorMeridianMm',0.5,@isnumeric);

% parse
p.parse(varargin{:})


%% Perform the calculation

% Determine polarAngle value for the vector that arises from the visual
% axis and intersects the optical axis
angleVisualToOpticalAxis = rad2deg(tan( p.Results.displacementOpticalAxisFromVisualAlongSuperiorMeridianMm / ...
    p.Results.displacementOpticalAxisFromVisualAlongNasalMeridianMm));

% Determine the length of the vector that arises from the visual axis and
% intersects the optical axis
distanceMmRetinaVisualToOpticalAxis = sqrt(p.Results.displacementOpticalAxisFromVisualAlongSuperiorMeridianMm^2 + ...
    p.Results.displacementOpticalAxisFromVisualAlongNasalMeridianMm^2);


end

