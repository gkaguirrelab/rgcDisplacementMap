function [ angleVisualToOpticalAxis, distanceDegVisualFieldVisualToOpticalAxis ] = visualFieldVectorVisualToOpticalAxis( varargin )
% Vector (angle and distance) relating visual -> optical axis of the eye
%
% Description:
%   Returns the polar angle and distance in degrees within the visual field
%   of the vector that arises at the visual axis of the eye and terminates
%   at the optical axis of the eye.
%
%   The coordinate system for polar angle is in degrees, and corresponds to
%       0=nasal; 90=superior; 180=temporal; 270=inferior
%   on the RETINAL field.
%
%   Note again the combination of units in the output:
%       angle - the polar angle ON THE RETINA of the line that connects the
%               visual to the optical axis of the eye
%       distance - degrees of visual field of the length of the line
%
%   The returned values are entirely dependent upon the settings of the two
%   optional analysis parameters. The default values are from Watson
%	2014 (in turn from Charman 1991, quoting Emsley 1952), and state that
%	the optical axis intersects the retina 1.5 mm nasal and 0.5 mm superior
%	to the visual axis.
%
% Inputs:
%   none
%
% Optional key/value pairs:
%   displacementOpticalAxisFromVisualAlongNasalMeridianMm
%   displacementOpticalAxisFromVisualAlongSuperiorMeridianMm
%
% Outputs:
%   angleVisualToOpticalAxis - Polar angle on the retinal field for the
%                           vector that connects the visual to the optical
%                           axis of the eye
%   distanceDegVisualFieldVisualToOpticalAxis - Distance (in degrees of
%                           visual angle of this vector).
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
% intersects the optical axis in units of degrees of visual field
distanceMmRetinaVisualToOpticalAxis = sqrt(p.Results.displacementOpticalAxisFromVisualAlongSuperiorMeridianMm^2 + ...
    p.Results.displacementOpticalAxisFromVisualAlongNasalMeridianMm^2);
distanceDegVisualFieldVisualToOpticalAxis = convert_mmRetina_to_degVisual(distanceMmRetinaVisualToOpticalAxis, angleVisualToOpticalAxis);


end

