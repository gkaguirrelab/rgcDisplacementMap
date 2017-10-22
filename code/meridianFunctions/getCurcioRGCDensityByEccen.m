function [rgcDensitySqDeg, supportPosDeg] = getCurcioRGCDensityByEccen(polarAngle)
% getCurcioRGCDensityByEccen(polarAngle)
%
% This routine returns the RGC density data reported in:
%
%   Curcio, Christine A., and Kimberly A. Allen. "Topography of ganglion
%   cells in human retina." Journal of comparative Neurology 300.1 (1990):
%   5-25.
% 
% Curcio and Allen obtained measurements of the density of all RGC classes
% within 6 human retinas at a set of positions relative to the fovea. These
% data were provided online in 2013.

% Here, we load these data and convert from mm to degrees.
%
% Inputs:
%   polarAngle - The desired angle of the density function on the retinal field.
%                (0=nasal;90=superior;180=temporal;270=inferior)
% Outputs:
%   rgcDensitySqDeg - the density (counts per square
%       degree) of RGCs at each of the positions
%   supportPosDeg - the positions (in degrees of visual angle) from the
%       fovea at which the cone density is defined
%

% Check the input
if sum([0 90 180 270]==polarAngle) ~= 1
    error('The Curcio RGC data are defined only for the cardinal meridians');
end

%% Load the RGC Density Data from Curcio et al 1990:
% Curcio and Allen obtained measurements of the density of RGCs
% within 6 human retinas at a set of positions relative to the fovea.
% Loading this matlab brings the variable "curcioRGCDensityPerSqMm" into
% memory
curcioRGCDataFile = ...
    fullfile([getpref('rgcDisplacementMap','LocalDataPath') , '/Curcio_1990_JCompNeurol_GanglionCellTopography/curcioRGCDensityPerSqMm.mat']);
load(curcioRGCDataFile);

meridianNames = {'nasal','superior','temporal','inferior'};

% Convert mm to deg, and mm^2 to deg^2
curcioRGCDensityPerSqDeg.support = ...
    convert_mm_to_deg(curcioRGCDensityPerSqMm.support);

for mm = 1:length(meridianNames)
curcioRGCDensityPerSqDeg.(meridianNames{mm}) = ...
    convert_mmSq_to_degSq(curcioRGCDensityPerSqDeg.support, curcioRGCDensityPerSqMm.(meridianNames{mm}) );
end

supportPosDeg = curcioRGCDensityPerSqDeg.support;
switch polarAngle
    case 0
        rgcDensitySqDeg = curcioRGCDensityPerSqDeg.nasal;
    case 90
        rgcDensitySqDeg = curcioRGCDensityPerSqDeg.superior;
    case 180
        rgcDensitySqDeg = curcioRGCDensityPerSqDeg.temporal;
    case 270
        rgcDensitySqDeg = curcioRGCDensityPerSqDeg.inferior;
end

end % function


