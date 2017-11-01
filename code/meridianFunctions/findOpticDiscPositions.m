
function opticDiscIndices = findOpticDiscPositions(regularSupportPosDeg, polarAngle, varargin)

%   opticDiscParameters - The model of the optic disc is govened by:
%       [ m, e, major, minor, theta]
%           m - meridian in degrees
%           e - eccentricity of the center of the optic disc in deg retina
%           major - length of the major axis in degrees retina
%           minor - length of the minor axis in degrees retina
%           theta - orientation of the ellipse in meridian degrees
%       Default values taken from Drasdo & Fowler 1974.

%% Parse input and define variables
p = inputParser;

% required input
p.addRequired('regularSupportPosDeg',@isnumeric);
p.addRequired('polarAngle',@isnumeric);

% Optional analysis params
p.addParameter('opticDiscParameters',struct('m',0,'e',5,'major',2,'minor',1,'theta',90),@isstruct);

% parse
p.parse(regularSupportPosDeg, polarAngle, varargin{:})

% Determine the distance between the center of the optic disc and
% each point in regularSupportPosDeg
distanceToOpticDiscCenterDegRetina = sqrt(...
    opticDiscParameters.e^2 + ...
    regularSupportPosDeg.^2 - ...
    2.*opticDiscParameters.e.*opticDiscParameters.e.*cos(deg2rad(polarAngle-opticDiscParameters.m)));

% Determine the angle (theta) between the center of the optic disc and
% each point in regularSupportPosDeg, relative to the major axis of the
% optic disc
angleWithOpticDiscCenter = rad2deg(arcsin((sin(deg2rad(polarAngle-opticDiscParameters.m))./regularSupportPosDeg).*distanceToOpticDiscCenterDegRetina))-opticDiscParameters.theta;


% The equation of an ellipse in polar coordinate form
r = (opticDisc.major*opticDisc.minor) ./ sqrt( (opticDisc.minor.*cos(deg2rad(angleWithOpticDiscCenter)).^2 + (opticDisc.major.*sin(deg2rad(angleWithOpticDiscCenter)).^2  )));



end