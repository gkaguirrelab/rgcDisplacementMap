function [ ] = displayRetinalImage(image, clim, pixelsPerDegree, maxEccentricity, colorBarString)
% A standardized display format for images of retinal data 
%
% Description:
%	This routine displays images up-down reversed with respect to the usual
%   convention used by Matlab, but corresponding to the usual orientation
%   of presentation of images of the fundus of the retina. Additionally,
%   points in the image that have a value of nan are displayed with black.
%
%   The expectation of the routine is that a figure with active axes
%   already exists.
%
% Inputs:
%   image                 - The image to be displayed
%   clim                  - A two element vector providing the min and max
%                           values to be used to defined the color bar
%   pixelsPerDegree       - Resolution of the image
%   maxEccentricity       - The maximum eccentricity represented in the
%                           image
%   colorBarString        - The string to be placed aside the color bar
%
% Outputs:
%   none
%

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
c.Label.String=colorBarString;
xlabel('Position [retinal deg] nasal --> temporal');
ylabel('Position [retinal deg] inferior --> superior');
xTickValuesInPixels = xticks;
xticklabels(string((xTickValuesInPixels.*(1/pixelsPerDegree)) + (-1*maxEccentricity)));
yticklabels(string(fliplr( (xTickValuesInPixels.*(1/pixelsPerDegree)) + (-1*maxEccentricity) - (xTickValuesInPixels(1)/pixelsPerDegree) )));


end

