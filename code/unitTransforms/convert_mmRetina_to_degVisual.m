function supportPosDegVisual = convert_mmRetina_to_degVisual(supportPosMmRetina)
% convert_mmRetina_to_degVisual
%
%   Converts mm on the retina from the visual axis to visual angle in
%   degrees from the visual axis. It is based on Appendix 6 of Watson 2014,
%   which in turn is derived from Drasdo & Fowler 1974.
%
% Input: 
%   supportPosMm  = Retinal postion(s) in milimeters. Either scalar value or
%                     vector input accepted.
%
% Output: 
%   supportPosDeg = Retinal postion(s) in degrees
%
% Sample Call:
%   supportPosMm  = 0:2:20;
%   supportPosDeg = convert_deg_to_mm(supportPosMm);
%

supportPosDegVisual = 3.556.*supportPosMmRetina + 0.05593.*(supportPosMmRetina.^2) - 0.007358.*(supportPosMmRetina.^3) +0.0003027.*(supportPosMmRetina.^4);

end

