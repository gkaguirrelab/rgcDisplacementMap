function imageMap = convertPolarMapToImageMap(polarMap, imRdim)

maxMapDensity = max(polarMap(:));
imP=polarMap'./maxMapDensity;
imR = PolarToIm (imP, 0, 1, imRdim, imRdim);
imageMap = imrotate(imR .* maxMapDensity,-90);

% Set values beyond the radius of the provided meridians to nans
distanceMap = zeros(size(imageMap));
distanceMap(ceil(imRdim/2),ceil(imRdim/2))=1;
distanceMap=bwdist(distanceMap,'euclidean');
imageMap(distanceMap>(imRdim/2))=nan;

end