% basicDemo.m
% Demonstrate the creation of a default RGC displacement map

%% Create the default displacement model.
% The variable rgcDisplacementByMeridian contains an m x p matrix, where m
% is number of meridians and p is the number of modeled eccentricity
% positions. The values specify the displacement of the RGC soma away from
% the fovea in units of retinal degrees. The next two variables are the
% support for the m and p dimensions, respectively.
[rgcDisplacementByMeridian, meridianAngleSupport, regularSupportPosDegRetina] = ...
    createDisplacementModel('verbose',true);


%% Convert from polar to image coordinates

% Define the image sample base for transform from polar coords. We reduce
% the resolution from the default of 0.01 retinal degrees per pixel to 10
% per pixel to speed things up.
displacementMapPixelsPerRetinalDeg = 10;
maxModeledEccentricity = max(regularSupportPosDegRetina);
imRdim = ...
    (maxModeledEccentricity * displacementMapPixelsPerRetinalDeg * 2) -1;

% This step actually makes the image map. If the 'imRdim' key/value pair is
% omitted, the map is returned with the same spatial resolution as the
% variable rgcDisplacementByMeridian. 
rgcDisplacementMap = ...
    convertPolarMapToImageMap(rgcDisplacementByMeridian,'imRdim',imRdim);


%% Display the map
displayRetinalImage(...
    rgcDisplacementMap, ...                     % the image to be displayed
    [0 ceil(max(max(rgcDisplacementMap)))], ... % min max of the color bar range
    displacementMapPixelsPerRetinalDeg, ...     % resolution of the image map
    maxModeledEccentricity, ...                 % the maximum eccentricity of the model
    'Displacement in retinal degrees');         % label for the color bar

