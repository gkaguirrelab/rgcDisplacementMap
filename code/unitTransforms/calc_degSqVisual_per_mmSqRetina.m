function degSqVisualPerMmSqRetinaRelativeToVisualAxis = calc_degSqVisual_per_mmSqRetina(supportPosMmRetinaRelativeToVisualAxis, polarAngle)
% Factor that expresses the deg^2 in the visual field per mm^2 retina
%
% Description:
%   Returns the square degree of visual angle per mm^2 on the retina as a
%   function of position in mm on the retina relative to the visual axis of
%   the eye. It is based on Appendix 6 of Watson 2014, which in turn is
%   derived from Drasdo & Fowler 1974. Drasdo & Fowler's equation relates
%   position on the retina to position in the visual field relative to the
%   optic axis of the eye. Watson 2014 observed that a correction
%   (dependent on polar angle) is necessary to correct the conversion to be
%   relative to the visual axis (which is otherwise assumed throughout this
%   toolbox).  Watson implemented an approximation to this correction; we
%   implement the full geometric correction here.
%
% Inputs:
%   supportPosMmRetinaRelativeToVisualAxis - retinal postion(s) in
%                           milimeters relative to visual axis of the eye.
%                           Either scalar value or vector input accepted.
%   polarAngle            - The polar angle of the meridian for which the
%                           conversion should be calculated. An arbitrary
%                           value between 0-360 is accepted. (0=nasal;
%                           90=superior; 180=temporal; 270=inferior)
%
% Outputs:
%   degSqVisualPerMmSqRetinaRelativeToVisualAxis - the area in square
%                           degrees in the visual field corresponding to a
%                           square mm of retina at the locations
%                           corresponding to
%                           supportPosMmRetinaRelativeToVisualAxis
%


%% Parse input and define variables
p = inputParser;

% required input
p.addRequired('supportPosMmRetinaRelativeToVisualAxis',@isnumeric);
p.addRequired('polarAngle',@isnumeric);

% parse
p.parse(supportPosMmRetinaRelativeToVisualAxis,polarAngle)

%% Perform the calculation

% Convert retinal mm to visual field degrees
supportPosDegVisualRelativeToVisualAxis = ...
    convert_mmRetina_to_degVisual(supportPosMmRetinaRelativeToVisualAxis, polarAngle);

% Obtain the area conversion at these locations
mmSqRetinaPerDegSqVisualRelativeToVisualAxis = ...
    calc_mmSqRetina_per_degSqVisual(supportPosDegVisualRelativeToVisualAxis, polarAngle);

% Return the reciprocal
degSqVisualPerMmSqRetinaRelativeToVisualAxis = 1 ./ mmSqRetinaPerDegSqVisualRelativeToVisualAxis;

end % main function
