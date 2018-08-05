function midgetFraction = calcDrasdoMidgetFractionByVisualEccen(supportPosDegVisual,f0,rm)
% Calculate the fraction of midget receptive fields
%
% Description:
%   Returns the fraction of the total receptive fields at given visual
%   field position that are midget receptive fields. The equation is taken
%   from Drasdo et al 2007.
%
% Inputs:
%   supportPosDegVisual   - A 1 x p vector that identifies the eccentricity
%                           positions (in visual degrees relative to the 
%                           visual axis / fovea) at which to perform the
%                           calculation.
%   f0                    - The midget fraction assigned to the fovea
%   rm                    - Decay parameter controlling the midget fraction
%                           function
%
% Outputs:
%   midgetFraction        - A 1 x p vector that contains the fraction of
%                           midget receptive fields at each eccentricity
%                           location.
%

midgetFraction = f0.*(1+(supportPosDegVisual./rm)).^-1;

end