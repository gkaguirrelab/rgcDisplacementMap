function supportPosMmRetina = convert_degVisual_to_mmRetina(supportPosDegVisual)
% convert_degVisual_to_mmRetina
%
%   Converts visual angle in degrees from the visual axis to mm on the
%   retina from the visual axis. It is based on Appendix 6 of Watson 2014,
%   which in turn is derived from Drasdo & Fowler 1974.
%
% Input: 
%   supportPosDeg = Retinal postion(s) in degrees. Either scalar value or
%                   vector input accepted.
%
% Output: 
%   supportPosMm  = Retinal postion(s) in millimeters
%
% Sample Call:
%   supportPosDeg = 0:2:20;
%   supportPosMm  = convert_deg_to_mm(supportPosDeg);
%

supportPosMmRetina = 0.268.*supportPosDegVisual + 0.0003427.*(supportPosDegVisual.^2) - 8.3309e-6.*(supportPosDegVisual.^3);

end

