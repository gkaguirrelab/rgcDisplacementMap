function validation_IndividualSubjectDensityFunctions(varargin)
% validation_IndividualSubjectDensityFunctions
%


%% Parse input and define variables
p = inputParser;

% Optional anaysis params
p.addParameter('cardinalMeridianNames',{'nasal','superior','temporal','inferior'},@iscell);
p.addParameter('cardinalMeridianAngles',[0 90 180 270],@isnumeric);
p.addParameter('cardinalMeridianPlotColors',{'r','b','g','k'},@iscell);
p.addParameter('subjectName', '29986A', @ischar);

% Optional display and ouput params
p.addParameter('verbose',true,@islogical);
p.addParameter('savePlots',true,@islogical);
p.addParameter('pathToPlotOutputDir','~/Desktop/rgcDisplacementMapPlots',@ischar);

% parse
p.parse(varargin{:})


%% Clean up
close all

% Define the cone and RGC density files to work with

coneDensityDataFileName = fullfile([getpref('rgcDisplacementMap','LocalDataPath') , '/Curcio_1990_JCompNeurol_HumanPhotoreceptorTopography/curcioRawConeDensity_',p.Results.subjectName,'.mat']);
rgcDensityDataFileName = fullfile([getpref('rgcDisplacementMap','LocalDataPath') , '/Curcio_1990_JCompNeurol_GanglionCellTopography/curcioRawRGCDensity_',p.Results.subjectName,'.mat']);

[ displacementMapDeg, fitParams, meridianAngles, rgcDisplacementEachMeridian, mRGC_cumulativeEachMeridian, mRF_cumulativeEachMeridian, convergenceEccen] = ...
    makeDisplacementMap( 'coneDensityDataFileName', coneDensityDataFileName, 'rgcDensityDataFileName', rgcDensityDataFileName, 'verbose', true );


end % function