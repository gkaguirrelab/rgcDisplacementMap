function mmSqRetinaPerDegSqVisualRelativeToVisualAxis = calc_mmSqRetina_per_degSqVisual(supportPosDegVisualRelativeToVisualAxis, polarAngle)
% Factor that expresses the mm^2 retina per deg^2 in the visual field
%
% Description:
%   Returns the mm^2 on the retina per square degree of visual angle as a
%   function of position in degrees visual angle relative to the visual
%   axis of the eye. It is based on Appendix 6 of Watson 2014, which in
%   turn is derived from Drasdo & Fowler 1974. Drasdo & Fowler's equation
%   relates position on the retina to position in the visual field relative
%   to the optic axis of the eye. Watson 2014 observed that a correction
%   (dependent on polar angle) is necessary to correct the conversion to be
%   relative to the visual axis (which is otherwise assumed throughout this
%   toolbox). Watson implemented an approximation to this correction; we
%   implement the full geometric correction here.
%
% Inputs:
%   supportPosDegVisualRelativeToVisualAxis - Retinal postion(s) in
%                           milimeters relative to visual axis of the eye.
%                           Either scalar value or vector input accepted.
%   polarAngle            - The polar angle of the meridian for which the
%                           conversion should be calculated. An arbitrary
%                           value between 0-360 is accepted. (0=nasal;
%                           90=superior; 180=temporal; 270=inferior)
%
% Outputs:
%   mmSqRetinaPerDegSqVisualRelativeToVisualAxis - the area in mm^2 on the
%                           retina for a degree of visual angle at the
%                           locations corresponding to
%                           supportPosDegVisualRelativeToVisualAxis
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
[ angleVisualToOpticalAxis, distanceDegVisualFieldVisualToOpticalAxis ] = ...
    visualFieldVectorVisualToOpticalAxis();

% Convert distances from visual axis to distances from optical axis
supportPosDegVisualRelativeToOpticalAxis = ...
    sqrt(...
    (supportPosDegVisualRelativeToVisualAxis-(distanceDegVisualFieldVisualToOpticalAxis.*cos(deg2rad(polarAngle - angleVisualToOpticalAxis)))).^2 + ...
    (distanceDegVisualFieldVisualToOpticalAxis .* sin(deg2rad(polarAngle - angleVisualToOpticalAxis))).^2 ...
    );

% perform the computation
mmSqRetinaPerDegSqVisualRelativeToVisualAxis = ...
    drasdoAndFowlerDegVisualSqToMmRetinaSq(supportPosDegVisualRelativeToOpticalAxis);

end % main function

% Local function
function mmSqRetinaPerDegSqVisualRelativeToOpticalAxis = drasdoAndFowlerDegVisualSqToMmRetinaSq(supportPosDegVisualRelativeToOpticalAxis)
mmSqRetinaPerDegSqVisualRelativeToOpticalAxis = ...
    0.0752+5.846e-5*supportPosDegVisualRelativeToOpticalAxis-1.064e-5*supportPosDegVisualRelativeToOpticalAxis.^2+4.116e-8*supportPosDegVisualRelativeToOpticalAxis.^3;
end