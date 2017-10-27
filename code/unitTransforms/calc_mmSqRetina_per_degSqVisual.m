function mmSqRetinaPerDegSqVisual = calc_mmSqRetina_per_degSqVisual(supportPosDegVisual)
% calc_mmSqRetina_per_degSqVisual
%
%   Returns the mm^2 on the retina per square degree of visual angle as a
%   function of position in degrees visual angle. It is based on Appendix 6
%   of Watson 2014, which in turn is derived from Drasdo & Fowler 1974.
%
% Inputs:
%   supportPosDegVisual - sample position in degrees 
%
% OutPuts:
%   mmSqRetinaPerDegSqVisual - the area in mm^2 on the retina for a degree
%     of visual angle at the locations corresponding to the
%     supportPosDegVisual
%

% Calculate the alpha conversion factor. It varies by eccentricity position
% Units of alpha = mm^2/degVisual^2
mmSqRetinaPerDegSqVisual = 0.0752+5.846e-5*supportPosDegVisual-1.064e-5*supportPosDegVisual.^2+4.116e-8*supportPosDegVisual.^3;

end

