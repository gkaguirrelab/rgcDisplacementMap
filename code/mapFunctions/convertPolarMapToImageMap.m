function imageMap = convertPolarMapToImageMap(polarMap, imRdim)

maxMapDensity = max(polarMap(:));
imP=polarMap'./maxMapDensity;
imR = PolarToIm (imP, 0, 1, imRdim, imRdim);
imageMap = imrotate(imR .* maxMapDensity,-90);

end