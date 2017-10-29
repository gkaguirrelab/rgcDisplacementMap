function supportPosDegVisualFieldRelativeToVisualAxis = convert_mmRetina_to_degVisual(supportPosMmRetinaRelativeToVisualAxis, polarAngle)
% convert_mmRetina_to_degVisual
%
%   Converts mm on the retina from the visual axis to visual angle in
%   degrees from the visual axis. It is based on Appendix 6 of Watson 2014.
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
%   supportPosMmRetinaRelativeToVisualAxis - retinal postion(s) in
%       milimeters relative to visual axis of the eye. Either scalar value
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
p.addRequired('supportPosMmRetinaRelativeToVisualAxis',@isnumeric);
p.addRequired('polarAngle',@isnumeric);

% parse
p.parse(supportPosMmRetinaRelativeToVisualAxis,polarAngle)


%% Adjust for displacement of the visual from the optical axis
% check the input
if polarAngle>360 || polarAngle<0
    error('Provide polarAngle between 0 and 360 degrees');
end

% Get the vector that connects the visual to optical axis on the retinal
% field
[ angleVisualToOpticalAxis, distanceMmRetinaVisualToOpticalAxis ] = retinalVectorVisualToOpticalAxis();

% Convert distances from visual axis to distances from optical axis
supportPosMmRetinaRelativeToOpticalAxis = ...
    sqrt(...
    (supportPosMmRetinaRelativeToVisualAxis-(distanceMmRetinaVisualToOpticalAxis.*cos(deg2rad(polarAngle - angleVisualToOpticalAxis)))).^2 + ...
    (distanceMmRetinaVisualToOpticalAxis .* sin(deg2rad(polarAngle - angleVisualToOpticalAxis))).^2 ...
    );

%% Perform the calculation
supportPosDegVisualFieldRelativeToVisualAxis = ...
    drasdoAndFowlerConversionRetinaToField(supportPosMmRetinaRelativeToVisualAxis - supportPosMmRetinaRelativeToOpticalAxis) + ...
    drasdoAndFowlerConversionRetinaToField(supportPosMmRetinaRelativeToOpticalAxis);

end % main function

% Local Drasdo & Fowler equation
function degreesVisualRelativeToOpticAxis = drasdoAndFowlerConversionRetinaToField(mmRetinaRelativeToOpticAxis)
degreesVisualRelativeToOpticAxis = 3.556.*mmRetinaRelativeToOpticAxis + 0.05593.*(mmRetinaRelativeToOpticAxis.^2) - 0.007358.*(mmRetinaRelativeToOpticAxis.^3) +0.0003027.*(mmRetinaRelativeToOpticAxis.^4);
end

