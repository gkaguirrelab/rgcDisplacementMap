function supportPosMmRetinaRelativeToVisualAxis = convert_degVisual_to_mmRetina(supportPosDegVisualRelativeToVisualAxis, polarAngle, varargin)
% convert_degVisual_to_mmRetina
%
%   Converts visual angle in degrees from the visual axis to mm on the
%   retina from the visual axis. It is based on Appendix 6 of Watson 2014.
%   The core of the conversion is from Drasdo & Fowler 1974. Drasdo &
%   Fowler's equation relates position on the retina to position in the
%   visual field relative to the optic axis of the eye. Watson 2014
%   observed that a correction (dependent on polar angle) is necessary to
%   correct the conversion to be relative to the visual axis (which is
%   otherwise assumed throughout this toolbox). Watson offered
%   approximation to the correction; here we implement the full geometric
%   correction.
%
%
% Input:
%   supportPosDegVisualRelativeToVisualAxis - visual field position in
%       degrees relative to visual axis of the eye. Either scalar value
%       or vector input accepted.
%   polarAngle - the polarAngle of the meridian for which the conversion
%       should be calculated. An arbitrary value between 0-360 is accepted.
%       (0=nasal;90=superior;180=temporal;270=inferior)
%
% Output:
%   supportPosDegVisualFieldRelativeToVisualAxis
%
% Optional:
%   displacementOpticalAxisFromVisualAlongNasalMeridianMm
%   displacementOpticalAxisFromVisualAlongSuperiorMeridianMm - From Watson
%       2014 (in turn from Charman 1991, quoting Emsley 1952), the optical
%       axis intersects the retina 1.5 mm nasal and 0.5 mm superior to the
%       visual axis.
%


%% Parse input and define variables
p = inputParser;

% required input
p.addRequired('supportPosDegVisualRelativeToVisualAxis',@isnumeric);
p.addRequired('polarAngle',@isnumeric);

% Optional analysis params
p.addParameter('displacementOpticalAxisFromVisualAlongNasalMeridianMm',1.5,@isnumeric);
p.addParameter('displacementOpticalAxisFromVisualAlongSuperiorMeridianMm',0.5,@isnumeric);

% parse
p.parse(supportPosDegVisualRelativeToVisualAxis,polarAngle,varargin{:})


%% Adjust for displacement of the visual from the optical axis
% check the input
if polarAngle>360 || polarAngle<0
    error('Provide polarAngle between 0 and 360 degrees');
end

% Determine polarAngle value for the vector that arises from the visual
% axis and intersects the optical axis
angleVisualToOpticalAxis = rad2deg(tan( p.Results.displacementOpticalAxisFromVisualAlongSuperiorMeridianMm / ...
    p.Results.displacementOpticalAxisFromVisualAlongNasalMeridianMm));

% Determine the length of the vector that arises from the visual axis and
% intersects the optical axis in units of degrees of visual field
distanceMmRetinaVisualToOpticalAxis = sqrt(p.Results.displacementOpticalAxisFromVisualAlongSuperiorMeridianMm^2 + ...
    p.Results.displacementOpticalAxisFromVisualAlongNasalMeridianMm^2);
distanceDegVisualFieldVisualToOpticalAxis = convert_mmRetina_to_degVisual(distanceMmRetinaVisualToOpticalAxis, angleVisualToOpticalAxis);

% Convert distances from visual axis to distances from optical axis
supportPosDegVisualRelativeToOpticalAxis = ...
    sqrt(...
    (supportPosDegVisualRelativeToVisualAxis-(distanceDegVisualFieldVisualToOpticalAxis.*cos(deg2rad(polarAngle - angleVisualToOpticalAxis)))).^2 + ...
    (distanceDegVisualFieldVisualToOpticalAxis .* sin(deg2rad(polarAngle - angleVisualToOpticalAxis))).^2 ...
    );

%% Perform the calculation
supportPosMmRetinaRelativeToVisualAxis = ...
    drasdoAndFowlerConversionFieldToRetina(supportPosDegVisualRelativeToVisualAxis - supportPosDegVisualRelativeToOpticalAxis) + ...
    drasdoAndFowlerConversionFieldToRetina(supportPosDegVisualRelativeToOpticalAxis);

end % main function

% Local Drasdo & Fowler equation
function mmRetinaRelativeToOpticAxis = drasdoAndFowlerConversionFieldToRetina(degreesVisualRelativeToOpticAxis)
mmRetinaRelativeToOpticAxis = 0.268.*degreesVisualRelativeToOpticAxis + 0.0003427.*(degreesVisualRelativeToOpticAxis.^2) - 8.3309e-6.*(degreesVisualRelativeToOpticAxis.^3);
end

