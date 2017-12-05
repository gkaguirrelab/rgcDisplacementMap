function opticDiscIndices = findOpticDiscPositions(regularSupportPosDegRetina, polarAngle, varargin)
% Model the spatial position of the optic disc on the retinal field
%
% Description:
%   Returns a vector that indicates if a given location in the input vector
%   of positions is located within the optic disc.
%
%   The default values were derived from the following sources:
%
%   K Rohrschneider. Determination of the Location of the Fovea on the
%		Fundus. Invest. Ophthalmol. Vis. Sci. 2004;45(9):3257-3258
%
%       This paper provides mean (across subject) values for the distance
%       of the optic disc center from the fovea. From these values I
%       derived the polar angle (5.6° visual) and the eccentricity
%       (15.5724° visual). The eccentricity is then converted to retinal
%       degrees.
%
% 	Drasdo & Fowler 1974
%
%       This paper provides the length (in mm) of the major (1.8 mm) and
%       minor (1.6 mm) axes of the optic disc. These are converted to
%       retinal degrees.
%
%   We assume that major axis of the optic disc is vertical (90°), but we
%   are aware that this can vary, and that in myopia the disc is noted to
%   tilt towards the fovea.
%
% Inputs:
%   regularSupportPosDegRetina - A 1 x p vector that contains the
%                           eccentricity in retinal degrees at which to
%                           determine the presence of the optic disc
%   polarAngle            - A scalar that specifies the polar angle (in
%                           degrees) for the meridian to be evalauted
%
% Optional key/value pairs:
%   opticDiscParameters   - The model of the optic disc is govened by
%                           [ m, e, major, minor, theta]:
%                           m - meridian in degrees relative to the nasal
%                               meridian
%                           e - eccentricity of the center of the optic
%                               disc in deg retina
%                           major - length of the major axis in degrees
%                               retina
%                           minor - length of the minor axis in degrees 
%                               retina
%                           theta - orientation of the ellipse in meridian
%                               degrees
%
% Outputs:
%   opticDiscIndices      - A 1 x p vector (correspoinding to the input
%                           regularSupportPosDegRetina) that indicates (1
%                           or 0) if the optic disc is present at this
%                           location
%


%% Parse input and define variables
p = inputParser;

% required input
p.addRequired('regularSupportPosDeg',@isnumeric);
p.addRequired('polarAngle',@isnumeric);

% Optional analysis params
p.addParameter('opticDiscParameters', ...
    struct(...
    'm',5.6, ...
    'e',convert_mmRetina_to_degRetina(convert_degVisual_to_mmRetina(15.5724, 5.6)),...
    'major',convert_mmRetina_to_degRetina(1.8),...
    'minor',convert_mmRetina_to_degRetina(1.6),...
    'theta',90),...
    @isstruct);

% parse
p.parse(regularSupportPosDegRetina, polarAngle, varargin{:})

opticDiscParameters = p.Results.opticDiscParameters;

% Determine the distance between the center of the optic disc and
% each point in regularSupportPosDeg
distanceToOpticDiscCenterDegRetina = sqrt(...
    opticDiscParameters.e^2 + ...
    regularSupportPosDegRetina.^2 - ...
    2.*opticDiscParameters.e.*regularSupportPosDegRetina.*cos(deg2rad(polarAngle-opticDiscParameters.m)));

% Determine the angle (alpha) between the line that connects the optic disc
% to each point in regularSupportPosDeg and the line that connects the
% optic disc to the center of the visual axis (fovea). As arcsin does not
% return values larger than 90°, we need to perform the calculation two
% different ways depending upon the lengths of the sides of the triangles
% involved.

% This calculation is valid if the position in regularSupportPosDeg is less
% than opticDiscParameters.e
alpha_case1 = rad2deg(asin(sin(deg2rad(polarAngle-opticDiscParameters.m)).*regularSupportPosDegRetina./distanceToOpticDiscCenterDegRetina));

% This calculation is valid if the position in regularSupportPosDeg is
% greater than opticDiscParameters.e
alpha_case2 = 180 - (polarAngle-opticDiscParameters.m) - ...
    rad2deg(asin(sin(deg2rad(polarAngle-opticDiscParameters.m))./distanceToOpticDiscCenterDegRetina.*opticDiscParameters.e));

% store the correct alpha for each regularSupportPosDeg
alpha=nan(size(alpha_case1));
alpha(regularSupportPosDegRetina <= opticDiscParameters.e) = ...
    alpha_case1(regularSupportPosDegRetina <= opticDiscParameters.e);
alpha(regularSupportPosDegRetina > opticDiscParameters.e) = ...
    alpha_case2(regularSupportPosDegRetina > opticDiscParameters.e);

% Now derive the angle (beta) between the center of the optic disc and
% each point in regularSupportPosDeg, relative to the major axis of the
% optic disc.
beta = opticDiscParameters.theta + opticDiscParameters.m - alpha;


% For each support position, determine if distanceToOpticDiscCenterDegRetina
% is less than the radius of the optic disc for that angle.
opticDiscIndices = distanceToOpticDiscCenterDegRetina < ...
                    ( ...
                       ((opticDiscParameters.major/2)*(opticDiscParameters.minor/2)) ./ sqrt( ((opticDiscParameters.minor/2).^2.*cos(deg2rad(beta)).^2 + ...
                       ((opticDiscParameters.major/2).^2.*sin(deg2rad(beta)).^2  ))) ...
                    );

% Handle the case in which a regularSupportPosDeg value is exactly on the
% optic disc center
if polarAngle == opticDiscParameters.m
    exactCenterIdx = find(regularSupportPosDegRetina==opticDiscParameters.e);
    if ~isempty(exactCenterIdx)
        opticDiscIndices(exactCenterIdx)=1;
    end
end

end % main function