function [ ] = displayRetinalImage(image, clim, pixelsPerDegree, maxEccentricity, titleString)
% displayRetinalImage - display image maps
%
%  This routine displays images up-down reversed with respect to the usual
%  convention used by Matlab, but corresponding to the usual orientation of
%  presentation of images of the fundus of the retina. Additionally, points
%  in the image that have a value of nan are displayed with black.


[nr,nc] = size(image);
pcolor([flipud(image) nan(nr,1); nan(1,nc+1)]);
caxis(clim);
shading flat;
set(gca, 'ydir', 'reverse');
axis square
set(gca,'TickLength',[0 0])
% Set the axis backgroud to dark gray
set(gcf,'Color',[1 1 1]); set(gca,'Color',[.75 .75 .75]); set(gcf,'InvertHardCopy','off');
c = colorbar;
c.Label.String=titleString;
xlabel('Position [retinal deg] nasal --> temporal');
ylabel('Position [retinal deg] inferior --> superior');
k=(xticks.*(1/pixelsPerDegree)) + (-1*maxEccentricity);
xticklabels(string(k))
yticklabels(string(fliplr(k)))


end

