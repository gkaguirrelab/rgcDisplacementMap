

%% Create the default displacement model.
% The variable rgcDisplacementByMeridian contains an m x p matrix, where m
% is number of meridians and p is the number of eccentricity positions
% (starting from zero) modeled. The values are the displacement in retinal
% degrees of the RGC soma away from the fovea. The next two variables are
% the support for the m and p dimensions, respectively.
[rgcDisplacementByMeridian, meridianAngleSupport, regularSupportPosDegRetina] = createDisplacementModel();


%% Display the displacement image

% Define the image sample base for transform from polar coords
displacementMapPixelsPerRetinalDeg = 10;	% resolution of the image map;
maxModeledEccentricity = max(regularSupportPosDegRetina);
imRdim = (maxModeledEccentricity * displacementMapPixelsPerRetinalDeg * 2) -1;

% Convert from polar to image coordinates
rgcDisplacementMap = convertPolarMapToImageMap(rgcDisplacementByMeridian, imRdim);

% Display the map
displayRetinalImage(...
    rgcDisplacementMap, ...                     % the image to be displayed
    [0,4], ...                                  % color bar range
    displacementMapPixelsPerRetinalDeg, ...     % resolution of the image map
    maxModeledEccentricity, ...                 % the maximum eccentricity of the model
    'Displaement in retinal degrees');                               % color bar

