function validation_WatsonFigs1and2(varargin)
% validation_WatsonFigs1and2
%
% This routine attempts to re-create Watson's Figures 1 and 2 starting
% From Curcio's raw data.

%% Parse input and define variables
p = inputParser;

% Optional anaysis params
p.addParameter('cardinalMeridianNames',{'nasal','superior','temporal','inferior'},@iscell);
p.addParameter('cardinalMeridianAngles',[0 90 180 270],@isnumeric);
p.addParameter('cardinalMeridianPlotColors',{'r','b','g','k'},@iscell);

% Optional display and ouput params
p.addParameter('verbose',true,@islogical);
p.addParameter('savePlots',true,@islogical);
p.addParameter('pathToPlotOutputDir','~/Desktop/rgcDisplacementMapPlots',@ischar);

% parse
p.parse(varargin{:})


%% Clean up
close all

%% Cone density functions

% make a figure
figure;
% loop over meridians
for mm = 1:length(p.Results.cardinalMeridianAngles)
    % Load the Curcio average data reported in the paper
    reportedConeDensityFile = fullfile([getpref('rgcDisplacementMap','LocalDataPath') , '/Curcio_1990_JCompNeurol_HumanPhotoreceptorTopography/curcioRawConeDensity_reportedAverage.mat']);
    [reportedConeDensitySqDegRetina, reportedConeDensitySupportDegRetina] = loadRawConeDensityByEccen(p.Results.cardinalMeridianAngles(mm), 'densityDataFileName', reportedConeDensityFile);

    % Convert to units of visual degrees
    reportedConeDensitySupportMmRetina = convert_degRetina_to_mmRetina(reportedConeDensitySupportDegRetina);
    reportedConeDensitySupportDegVisual = convert_mmRetina_to_degVisual(reportedConeDensitySupportMmRetina, p.Results.cardinalMeridianAngles(mm));
    reportedConeDensitySqMmRetina = reportedConeDensitySqDegRetina .* calc_degSqRetina_per_mmSqRetina();
    reportedConeDensitySqDegVisual = reportedConeDensitySqMmRetina ./ calc_degSqVisual_per_mmSqRetina(reportedConeDensitySupportMmRetina, p.Results.cardinalMeridianAngles(mm));

    % Plot the mRF:cone ratio    
    [ midgetRFDensitySqDegVisual ] = calcWatsonMidgetRFDensityByEccenDegVisual(reportedConeDensitySupportDegVisual, p.Results.cardinalMeridianAngles(mm));
%    loglog(reportedConeDensitySupportDegVisual, midgetRFDensitySqDegVisual ./ reportedConeDensitySqDegVisual, '-', 'Color', p.Results.cardinalMeridianPlotColors{mm});
    
    loglog(reportedConeDensitySupportDegVisual, reportedConeDensitySqDegVisual, '-', 'Color', p.Results.cardinalMeridianPlotColors{mm});
    hold on
    ylim([100 20000]);
    xlabel('degrees visual field');
    ylabel('counts/visual deg2');
    title('Watson 2014 Figure 1');    
end % loop over meridians
legend(p.Results.cardinalMeridianNames);


%% RGC density functions

% make a figure
figure;
% loop over meridians
for mm = 1:length(p.Results.cardinalMeridianAngles)
    % Load the Curcio average data reported in the paper
    reportedRGCDensityFile = fullfile([getpref('rgcDisplacementMap','LocalDataPath') , '/Curcio_1990_JCompNeurol_GanglionCellTopography/curcioRawRGCDensity_reportedAverage.mat']);
    [reportedRGCDensitySqDegRetina, reportedRGCDensitySupportDegRetina] = loadRawRGCDensityByEccen(p.Results.cardinalMeridianAngles(mm), 'densityDataFileName', reportedRGCDensityFile);

    % Convert to units of visual degrees
    reportedRGCDensitySupportMmRetina = convert_degRetina_to_mmRetina(reportedRGCDensitySupportDegRetina);
    reportedRGCDensitySupportDegVisual = convert_mmRetina_to_degVisual(reportedRGCDensitySupportMmRetina, p.Results.cardinalMeridianAngles(mm));
    reportedRGCDensitySqMmRetina = reportedRGCDensitySqDegRetina .* calc_degSqRetina_per_mmSqRetina();
    reportedRGCDensitySqDegVisual = reportedRGCDensitySqMmRetina ./ calc_degSqVisual_per_mmSqRetina(reportedRGCDensitySupportMmRetina, p.Results.cardinalMeridianAngles(mm));
    
    loglog(reportedRGCDensitySupportDegVisual, reportedRGCDensitySqDegVisual, '-', 'Color', p.Results.cardinalMeridianPlotColors{mm});
    hold on
        ylim([1 3000]);
    xlabel('degrees visual field');
    ylabel('counts/visual deg2');
    title('Watson 2014 Figure 2');    
end % loop over meridians
legend(p.Results.cardinalMeridianNames);



end % function