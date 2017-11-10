function supportPosMmRetinaRelativeToVisualAxis = convert_degVisual_to_mmRetina(supportPosDegVisualRelativeToVisualAxis, polarAngle)
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


%% Parse input and define variables
p = inputParser;

% required input
p.addRequired('supportPosDegVisualRelativeToVisualAxis',@isnumeric);
p.addRequired('polarAngle',@isnumeric);

% parse
p.parse(supportPosDegVisualRelativeToVisualAxis,polarAngle)


%% Adjust for displacement of the visual from the optical axis
% check the input
if polarAngle>360 || polarAngle<0
    error('Provide polarAngle between 0 and 360 degrees');
end

% Get the vector that connects the visual to optical axis
[ angleVisualToOpticalAxis, distanceDegVisualFieldVisualToOpticalAxis ] = visualFieldVectorVisualToOpticalAxis();

% Convert distances from visual axis to distances from optical axis
supportPosDegVisualRelativeToOpticalAxis = ...
    sqrt(...
    (supportPosDegVisualRelativeToVisualAxis-(distanceDegVisualFieldVisualToOpticalAxis.*cos(deg2rad(polarAngle - angleVisualToOpticalAxis)))).^2 + ...
    (distanceDegVisualFieldVisualToOpticalAxis .* sin(deg2rad(polarAngle - angleVisualToOpticalAxis))).^2 ...
    );

%% Perform the calculation
supportPosMmRetinaRelativeToVisualAxis = ...
    ( sign(supportPosDegVisualRelativeToVisualAxis - supportPosDegVisualRelativeToOpticalAxis) .* ...
    drasdoAndFowlerConversionFieldToRetina(abs(supportPosDegVisualRelativeToVisualAxis - supportPosDegVisualRelativeToOpticalAxis)) ) + ...
    drasdoAndFowlerConversionFieldToRetina(supportPosDegVisualRelativeToOpticalAxis);


end % main function


% Local Drasdo & Fowler equations

function mmRetinaRelativeToOpticAxis = drasdoAndFowlerConversionFieldToRetina(degreesVisualRelativeToOpticAxis)
%
% In Appendix 6, Watson (2014) provides polynomial approximations to the
% plots created by Drasdo & Fowler for the forward and inverse
% transformation of mm on the retina to position in the visual field. We
% find that these approximations result in discrepancies on the order of a
% tenth of a millimeter when we attempt to recover mm after passing through
% the forward and inverse transform. Consequently, we take the Watson 2014
% approximation for the forward transform (mm retina --> deg visual field),
% but for the inverse transform (deg visual field --> mm retina) we perform
% a spline fit to the inverted output of the forward transform. This
% provides a result that is closer to fully invertible.
%
supportPosMmRetinaRelativeToOpticalAxis=0:0.01:25;
calculatedPosDegFieldRelativeToOpticalAxis = ...
    drasdoAndFowlerConversionRetinaToField(supportPosMmRetinaRelativeToOpticalAxis);
mmRetinaRelativeToOpticAxis = ...
    spline(calculatedPosDegFieldRelativeToOpticalAxis,supportPosMmRetinaRelativeToOpticalAxis,degreesVisualRelativeToOpticAxis);
end

function degreesVisualRelativeToOpticAxis = drasdoAndFowlerConversionRetinaToField(mmRetinaRelativeToOpticAxis)
degreesVisualRelativeToOpticAxis = 3.556.*mmRetinaRelativeToOpticAxis + 0.05593.*(mmRetinaRelativeToOpticAxis.^2) - 0.007358.*(mmRetinaRelativeToOpticAxis.^3) +0.0003027.*(mmRetinaRelativeToOpticAxis.^4);
end
