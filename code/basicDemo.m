

%% Create the default displacement model.
% The variable rgcDisplacementByMeridian contains an m x p matrix, where m
% is number of meridians and p is the number of modeled eccentricity
% positions. The values are the displacement in retinal
% degrees of the RGC soma away from the fovea. The next two variables are
% the support for the m and p dimensions, respectively.
[rgcDisplacementByMeridian, meridianAngleSupport, regularSupportPosDegRetina] = ...
    createDisplacementModel('verbose',true);


%% Convert from polar to image coordinates
% Define the image sample base for transform from polar coords. We reduce
% the resolution a bit to speed things up.
displacementMapPixelsPerRetinalDeg = 10;
maxModeledEccentricity = max(regularSupportPosDegRetina);
imRdim = (maxModeledEccentricity * displacementMapPixelsPerRetinalDeg * 2) -1;
rgcDisplacementMap = ...
    convertPolarMapToImageMap(rgcDisplacementByMeridian,'imRdim',imRdim);


%% Display the map
displayRetinalImage(...
    rgcDisplacementMap, ...                     % the image to be displayed
    [0 4], ...                                  % min max of the color bar range
    displacementMapPixelsPerRetinalDeg, ...     % resolution of the image map
    maxModeledEccentricity, ...                 % the maximum eccentricity of the model
    'Displacement in retinal degrees');         % label for the color bar

