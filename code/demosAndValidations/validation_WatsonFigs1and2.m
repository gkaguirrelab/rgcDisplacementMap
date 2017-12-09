function validation_WatsonFigs1and2(varargin)
% Attempt to re-create Figures 1 and 2 from Watson 2014 JoV
%
% Description:
%   Figures 1 and 2 of Watson 2014 JoV present cone and RGC densities as a
%   function of eccentricity and cardinal meridian, with the values derived
%   from the data for the Curcio 1990 papers. Here, we attempt to re-create
%   Watson's figures starting from the raw data from Curcio. This tests if
%   we are implementing unit conversions in a manner that matches Watson.
%

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
    [reportedConeDensitySqDegRetina, reportedConeDensitySupportDegRetina] = loadRawConeDensityByEccen(p.Results.cardinalMeridianAngles(mm), 'coneDensityDataFileName', reportedConeDensityFile);

    % Convert to units of visual degrees
    reportedConeDensitySupportMmRetina = convert_degRetina_to_mmRetina(reportedConeDensitySupportDegRetina);
    reportedConeDensitySupportDegVisual = convert_mmRetina_to_degVisual(reportedConeDensitySupportMmRetina, p.Results.cardinalMeridianAngles(mm));
    reportedConeDensitySqMmRetina = reportedConeDensitySqDegRetina .* calc_degSqRetina_per_mmSqRetina();
    reportedConeDensitySqDegVisual = reportedConeDensitySqMmRetina ./ calc_degSqVisual_per_mmSqRetina(reportedConeDensitySupportMmRetina, p.Results.cardinalMeridianAngles(mm));

    % Plot the density functions
    loglog(reportedConeDensitySupportDegVisual, reportedConeDensitySqDegVisual, '-', 'Color', p.Results.cardinalMeridianPlotColors{mm});
    hold on
    ylim([100 20000]);
    xlim([0.1 100]);
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
    [reportedRGCDensitySqDegRetina, reportedRGCDensitySupportDegRetina] = loadRawRGCDensityByEccen(p.Results.cardinalMeridianAngles(mm), 'rgcDensityDataFileName', reportedRGCDensityFile);

    % Convert to units of visual degrees
    reportedRGCDensitySupportMmRetina = convert_degRetina_to_mmRetina(reportedRGCDensitySupportDegRetina);
    reportedRGCDensitySupportDegVisual = convert_mmRetina_to_degVisual(reportedRGCDensitySupportMmRetina, p.Results.cardinalMeridianAngles(mm));
    reportedRGCDensitySqMmRetina = reportedRGCDensitySqDegRetina .* calc_degSqRetina_per_mmSqRetina();
    reportedRGCDensitySqDegVisual = reportedRGCDensitySqMmRetina ./ calc_degSqVisual_per_mmSqRetina(reportedRGCDensitySupportMmRetina, p.Results.cardinalMeridianAngles(mm));
    
    % Plot the density functions
    loglog(reportedRGCDensitySupportDegVisual, reportedRGCDensitySqDegVisual, '-', 'Color', p.Results.cardinalMeridianPlotColors{mm});
    hold on
    ylim([1 3000]);
    xlim([0.1 100]);
    xlabel('degrees visual field');
    ylabel('counts/visual deg2');
    title('Watson 2014 Figure 2');    
end % loop over meridians
legend(p.Results.cardinalMeridianNames);



end % function