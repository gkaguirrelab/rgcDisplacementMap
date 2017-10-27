function validation_CurcioAverages(varargin)
% validation_CurcioAverages
%
% This routine compares the average cone and RGC densities reported in
% Christine Curcio's two 1990 J Comp Neurology papers to the averages that
% we compute from the raw individual data that she has provided.
%
% We find that the reported RGC averages in the Curcio paper are not equal
% to the average of the individual subject values. It seems likely
% that the average values reported in the paper were derived from an
% intermediate modeling step that we cannot currently reconstruct.

%% Parse input and define variables
p = inputParser;

% Optional anaysis params
p.addParameter('cardinalMeridianNames',{'nasal','superior','temporal','inferior'},@iscell);
p.addParameter('cardinalMeridianAngles',[0 90 180 270],@isnumeric);

% Optional display and ouput params
p.addParameter('verbose',true,@islogical);
p.addParameter('savePlots',true,@islogical);
p.addParameter('pathToPlotOutputDir','~/Desktop/rgcDisplacementMapPlots',@ischar);

% parse
p.parse(varargin{:})


%% Clean up
close all

%% RGC density functions


% make a figure
figure;

% loop over meridians
for mm = 1:length(p.Results.cardinalMeridianAngles)
    calculateRGCdensityFile = fullfile([getpref('rgcDisplacementMap','LocalDataPath') , '/Curcio_1990_JCompNeurol_GanglionCellTopography/curcioRawRGCDensity_computedAverage.mat']);
    reportedRGCdensityFile = fullfile([getpref('rgcDisplacementMap','LocalDataPath') , '/Curcio_1990_JCompNeurol_GanglionCellTopography/curcioRawRGCDensity_reportedAverage.mat']);
    
    [calculatedRGCDensitySqDegRetina, calculateRGCdensitySupportDegRetina] = loadRawRGCDensityByEccen(p.Results.cardinalMeridianAngles(mm), 'densityDataFileName', calculateRGCdensityFile);
    [reportedRGCDensitySqDegRetina, reportedRGCdensitySupportDegRetina] = loadRawRGCDensityByEccen(p.Results.cardinalMeridianAngles(mm), 'densityDataFileName', reportedRGCdensityFile);
    
    subplot(2,2,mm);
    plot(calculateRGCdensitySupportDegRetina, calculatedRGCDensitySqDegRetina, 'xb');
    hold on
    plot(reportedRGCdensitySupportDegRetina, reportedRGCDensitySqDegRetina, '.r');
    plot(calculateRGCdensitySupportDegRetina, calculatedRGCDensitySqDegRetina, '-b');
    xlabel('retinal degrees');
    ylabel('counts/mm2 retina');
    legend('computed','reported')
    title(p.Results.cardinalMeridianNames{mm});    
end % loop over meridians


%% Cone density functions

% make a figure
figure;

% loop over meridians
for mm = 1:length(p.Results.cardinalMeridianAngles)
    calculateRGCdensityFile = fullfile([getpref('rgcDisplacementMap','LocalDataPath') , '/Curcio_1990_JCompNeurol_HumanPhotoreceptorTopography/curcioRawConeDensity_computedAverage.mat']);
    reportedRGCdensityFile = fullfile([getpref('rgcDisplacementMap','LocalDataPath') , '/Curcio_1990_JCompNeurol_HumanPhotoreceptorTopography/curcioRawConeDensity_reportedAverage.mat']);
    
    [calculatedRGCDensitySqDegRetina, calculateRGCdensitySupportDegRetina] = loadRawRGCDensityByEccen(p.Results.cardinalMeridianAngles(mm), 'densityDataFileName', calculateRGCdensityFile);
    [reportedRGCDensitySqDegRetina, reportedRGCdensitySupportDegRetina] = loadRawRGCDensityByEccen(p.Results.cardinalMeridianAngles(mm), 'densityDataFileName', reportedRGCdensityFile);
    
    subplot(2,2,mm);
    plot(calculateRGCdensitySupportDegRetina, calculatedRGCDensitySqDegRetina, 'xb');
    hold on
    plot(reportedRGCdensitySupportDegRetina, reportedRGCDensitySqDegRetina, '.r');
    plot(calculateRGCdensitySupportDegRetina, calculatedRGCDensitySqDegRetina, '-b');
    xlabel('retinal degrees');
    ylabel('counts/mm2 retina');
    legend('computed','reported')
    title(p.Results.cardinalMeridianNames{mm});    
end % loop over meridians


end % function