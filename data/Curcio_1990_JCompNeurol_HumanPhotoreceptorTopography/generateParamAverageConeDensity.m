function generateParamAverageConeDensity(varargin)
% generateParamAverageConeDensity
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
p.addParameter('subjectNames',...
    {'E8746R', '9387L', '6884R', '59784L', '483848L', '44985L', '29986R' ,'29986L'}, ...
    @iscell);

% parse
p.parse(varargin{:})

outputFileNameStem = 'curcioRawConeDensity_';

subjectName = 'paramAverage';

%% Clean up
close all

curcioParamAverageConeDensity = struct;

for mm = 1:4 %length(p.Results.cardinalMeridianAngles)
    figure
    linearCoefficients=[];
    for ss = 1:length(p.Results.subjectNames)
        % Load the Curcio average data reported in the paper
        reportedConeDensityFile = fullfile([getpref('rgcDisplacementMap','LocalDataPath') , '/Curcio_1990_JCompNeurol_HumanPhotoreceptorTopography/curcioRawConeDensity_',p.Results.subjectNames{ss},'.mat']);
        [reportedConeDensitySqDegRetina, reportedConeDensitySupportDegRetinaInitial] = loadRawConeDensityByEccen(p.Results.cardinalMeridianAngles(mm), 'densityDataFileName', reportedConeDensityFile);

        % Give the zero eccentricity a slightly non-zero value to handle the log10
        reportedConeDensitySupportDegRetinaInitial(1) = 10.^(...
            log10(reportedConeDensitySupportDegRetinaInitial(2)) + ...
            log10(diff(reportedConeDensitySupportDegRetinaInitial(1:2))));

        % Remove nans
        nanIdx = isnan(reportedConeDensitySqDegRetina);
        if sum(nanIdx) ~=0
            reportedConeDensitySupportDegRetina = reportedConeDensitySupportDegRetinaInitial(~nanIdx);
            reportedConeDensitySqDegRetina = reportedConeDensitySqDegRetina(~nanIdx);
        end
        
        linearCoefficients(ss,:) = polyfit(log10(reportedConeDensitySupportDegRetina), log10(reportedConeDensitySqDegRetina), 1);
        yFit = polyval(linearCoefficients(ss,:), log10(reportedConeDensitySupportDegRetina));
        plot(log10(reportedConeDensitySupportDegRetina),log10(reportedConeDensitySqDegRetina),'.')
        hold on
        plot(log10(reportedConeDensitySupportDegRetina),yFit,'-')
    end % subjects
    
    medianLinearCoef = median(linearCoefficients);
    paramAverageConeDensitySqDegRetina = 10.^polyval(medianLinearCoef, log10(reportedConeDensitySupportDegRetinaInitial));
    reportedConeDensitySupportDegRetina(1)=0;
    curcioParamAverageConeDensity.support = reportedConeDensitySupportDegRetinaInitial;
    curcioParamAverageConeDensity.(lower(p.Results.cardinalMeridianNames{mm})) = paramAverageConeDensitySqDegRetina;

end % meridians

curcioParamAverageConeDensity.meta.subjectName = subjectName;
curcioParamAverageConeDensity.meta.densityUnits = 'counts/retinalDeg2';
curcioParamAverageConeDensity.meta.supportUnits = 'retinalDegree';
save([outputFileNameStem subjectName], 'curcioParamAverageConeDensity');


end % function

