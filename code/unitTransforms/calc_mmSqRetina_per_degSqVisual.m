function mmSqRetinaPerDegSqVisualRelativeToVisualAxis = calc_mmSqRetina_per_degSqVisual(supportPosDegVisualRelativeToVisualAxis, polarAngle)
% calc_mmSqRetina_per_degSqVisual
%
%   Returns the mm^2 on the retina per square degree of visual angle as a
%   function of position in degrees visual angle relative to the visual
%   axis of the eye. It is based on Appendix 6 of Watson 2014, which in
%   turn is derived from Drasdo & Fowler 1974. Drasdo &
%   Fowler's equation relates position on the retina to position in the
%   visual field relative to the optic axis of the eye. Watson 2014
%   observed that a correction (dependent on polar angle) is necessary to
%   correct the conversion to be relative to the visual axis (which is
%   otherwise assumed throughout this toolbox). It is not clear from
%   Watson's paper if he accounted for this in his implementation of this
%   equation; here we implement the full geometric correction.

% Inputs:
%   supportPosDegVisualRelativeToVisualAxis - retinal postion(s) in
%       milimeters relative to visual axis of the eye. Either scalar value
%       or vector input accepted.
%   polarAngle - the polarAngle of the meridian for which the conversion
%       should be calculated. An arbitrary value between 0-360 is accepted.
%       (0=nasal;90=superior;180=temporal;270=inferior)
%
% OutPuts:
%   mmSqRetinaPerDegSqVisualRelativeToVisualAxis - the area in mm^2 on the
%       retina for a degree of visual angle at the locations corresponding
%       to supportPosDegVisualRelativeToVisualAxis
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


% perform the computation
mmSqRetinaPerDegSqVisualRelativeToVisualAxis = ...
    drasdoAndFowlerAreaEquation(supportPosDegVisualRelativeToVisualAxis - supportPosDegVisualRelativeToOpticalAxis) + ...
    drasdoAndFowlerAreaEquation(supportPosDegVisualRelativeToOpticalAxis);

end % main function

% Local function
function mmSqRetinaPerDegSqVisualRelativeToOpticalAxis = drasdoAndFowlerAreaEquation(supportPosDegVisualRelativeToOpticalAxis)
mmSqRetinaPerDegSqVisualRelativeToOpticalAxis = ...
    0.0752+5.846e-5*supportPosDegVisualRelativeToOpticalAxis-1.064e-5*supportPosDegVisualRelativeToOpticalAxis.^2+4.116e-8*supportPosDegVisualRelativeToOpticalAxis.^3;
end