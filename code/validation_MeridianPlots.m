function validation_MeridianPlots(varargin)
% Create plots that illustrate the behavior of the model along the cardinal
% merdians


%% Parse input and define variables
p = inputParser; p.KeepUnmatched = true;

% Optional anaysis params
p.addParameter('sampleResolutionDegrees',0.01,@isnumeric);
p.addParameter('maxModeledEccentricity',30,@isnumeric);
p.addParameter('meridianAngleResolutionDeg',90,@isnumeric);
p.addParameter('displacementMapPixelsPerDeg',10,@isnumeric);
p.addParameter('pathToPlotOutputDir','~/Desktop/rgcDisplacementMapPlots',@ischar);

% Optional display and ouput params
p.addParameter('verbose',true,@islogical);
p.addParameter('savePlots',true,@islogical);

% parse
p.parse(varargin{:})


%% Setup
close all

% check or make a directory for output
if p.Results.savePlots
    if exist(p.Results.pathToPlotOutputDir,'dir')==0
        mkdir(p.Results.pathToPlotOutputDir);
    end
end

% Prepare the regular eccentricity support base
regularSupportPosDeg = ...
    0:p.Results.sampleResolutionDegrees:p.Results.maxModeledEccentricity;

% make the displacement map for the cardinal meridians
[ ~, fitParams, meridianAngles, rgcDisplacementEachMeridian, mRGC_cumulativeEachMeridian, mRF_cumulativeEachMeridian ] = ...
    makeDisplacementMap(...
    'sampleResolutionDegrees', p.Results.sampleResolutionDegrees, ...
    'maxModeledEccentricity', p.Results.maxModeledEccentricity, ...
    'meridianAngleResolutionDeg', p.Results.meridianAngleResolutionDeg, ...
    'displacementMapPixelsPerDeg', p.Results.displacementMapPixelsPerDeg, ...
    'verbose', p.Results.verbose);

% Plot the cumulatives and displacements across meridians
figHandle = figure();
for mm = 1:length(meridianAngles)
    
    % plot the displacement
    subplot(length(meridianAngles),2,mm*2);
    plot(regularSupportPosDeg,rgcDisplacementEachMeridian(mm,:),'-r')
    axis off;
    ylim([-.5 3.0]);
    if mm == length(meridianAngles)
        axis on;
        xlabel('eccentricity [deg]');
        ylabel('RGC displacement [deg]');
    end
    
    % Plot the cumulative functions
    subplot(length(meridianAngles),2,mm*2-1);
    plot(regularSupportPosDeg,mRGC_cumulativeEachMeridian(mm,:),'-k')
    axis off;
    if mm == length(meridianAngles)
        axis on;
        xlabel('eccentricity [deg]');
        ylabel('cells per sector');
    end
    hold on
    plot(regularSupportPosDeg,mRF_cumulativeEachMeridian(mm,:),'-b')
    ylim([0 8e5]);
    hold off
    drawnow
end
if p.Results.savePlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDir,'cumulativeAndDisplaceByMeridian.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end


% Show the cone --> mRF model
[ ~, figHandle ] = developMidgetRFFractionModel( 'makePlots', true );
if p.Results.savePlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDir,'cone --> mRF RatioModel.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end

% Show the RGC --> mRGC model
[ ~, figHandle ] = developMidgetRGCFractionModel( 'makePlots', true );
if p.Results.savePlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDir,'RGC --> mRGC RatioModel.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end

% Plot the spline cone density fits
figHandle = figure();
cardinalMeridianAngles=[0 90 180 270];
meridianColors={'r','b','g','k'};
subplot(1,2,1)
for mm=1:4
    [coneDensitySqDeg, coneNativeSupportPosDeg] = getCurcioConeDensityByEccen(cardinalMeridianAngles(mm));
    dispConeNativeSupportPosDeg=coneNativeSupportPosDeg;
    dispConeNativeSupportPosDeg(1)=1e-2;
    [coneDensityFit] = getSplineFitToConeDensity(cardinalMeridianAngles(mm));
    loglog(dispConeNativeSupportPosDeg,coneDensitySqDeg,'x','Color',meridianColors{mm});
    xlim([1e-2,70]);
    hold on
    loglog(dispConeNativeSupportPosDeg,coneDensityFit(coneNativeSupportPosDeg),'-','Color',meridianColors{mm});
    xlabel('log10 Eccentricity [deg]');
    ylabel('log10 Cone density [counts / deg2]');
end
legend(num2str(cardinalMeridianAngles),'Location','southwest')

subplot(1,2,2)
interpolarMeridianAngles=[45 135 225 315];
for mm=1:4
    [coneDensityFit] = getSplineFitToConeDensity(interpolarMeridianAngles(mm));
    loglog(dispConeNativeSupportPosDeg,coneDensityFit(coneNativeSupportPosDeg),'-','Color',meridianColors{mm});
    xlim([1e-2,70]);
    xlabel('log10 Eccentricity [deg]');
    ylabel('log10 Cone density [counts / deg2]');
    hold on
end
legend(num2str(interpolarMeridianAngles),'Location','southwest')
if p.Results.savePlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDir,'coneDensitySplineModel.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end


% Plot the spline RGC density fits
figHandle = figure();
cardinalMeridianAngles=[0 90 180 270];
meridianColors={'r','b','g','k'};
subplot(1,2,1)
for mm=1:4
    [RGCDensitySqDeg, RGCNativeSupportPosDeg] = getCurcioRGCDensityByEccen(cardinalMeridianAngles(mm));
    dispRGCNativeSupportPosDeg=RGCNativeSupportPosDeg;
    dispRGCNativeSupportPosDeg(1)=1e-2;
    [RGCDensityFit] = getSplineFitToRGCDensity(cardinalMeridianAngles(mm));
    loglog(dispRGCNativeSupportPosDeg,RGCDensitySqDeg,'x','Color',meridianColors{mm});
    xlim([1e-2,70]);
    hold on
    regularSupportPosDeg=1e-2:0.01:70;
    loglog(regularSupportPosDeg,RGCDensityFit(regularSupportPosDeg),'-','Color',meridianColors{mm});
    xlabel('log10 Eccentricity [deg]');
    ylabel('log10 RGC density [counts / deg2]');
end
legend(num2str(cardinalMeridianAngles),'Location','southwest')

subplot(1,2,2)
interpolarMeridianAngles=[45 135 225 315];
for mm=1:4
    [RGCDensityFit] = getSplineFitToRGCDensity(interpolarMeridianAngles(mm));
    loglog(regularSupportPosDeg,RGCDensityFit(regularSupportPosDeg),'-','Color',meridianColors{mm});
    xlim([1e-2,70]);
    xlabel('log10 Eccentricity [deg]');
    ylabel('log10 RGC density [counts / deg2]');
    hold on
end
legend(num2str(interpolarMeridianAngles),'Location','southwest')
if p.Results.savePlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDir,'rgcDensitySplineModel.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end


%% Plot the mRGC fraction for the cardinal meridians
figHandle = figure();
cardinalMeridianAngles=[0 90 180 270];
meridianColors={'r','b','g','k'};

for mm = 1:4
    
    meridianIdx = find(meridianAngles==cardinalMeridianAngles(mm),1);
    
    % Load the RGC Density Data from Curcio and Allen 1990:
    [ RGCDensitySqDeg, RGCNativeSupportPosDeg ] = getCurcioRGCDensityByEccen( cardinalMeridianAngles(mm) );
    % remove nan values
    isvalididx=find(~isnan(RGCDensitySqDeg)  );
    RGCNativeSupportPosDeg = RGCNativeSupportPosDeg(isvalididx);
    RGCDensitySqDeg = RGCDensitySqDeg(isvalididx);
    
    % Plot Watson's midget fraction
    subplot(1,2,1);
    f0 = 0.8928; rm = 41.03; % Watson's values
    midgetFraction_watson = calcWatsonMidgetFractionByEccen(RGCNativeSupportPosDeg,f0,rm);
    plot(RGCNativeSupportPosDeg,midgetFraction_watson,'-k');
    hold on
    xlabel('eccentricity deg');
    ylabel('midget fraction');
    ylim([0 1]);
    xlim([0 40]);
    title('Watson''s midget fraction (from Drasdo)');
    pbaspect([2 1 1]);
    
    % Plot our midget fraction
    subplot(1,2,2);
    [ ~, midgetFraction_ours ] = transformRGCToMidgetRGCDensity( RGCNativeSupportPosDeg, RGCDensitySqDeg, 'recipFitParams', fitParams(meridianIdx,3:5) );
    plot(RGCNativeSupportPosDeg,midgetFraction_ours,'-','Color',meridianColors{mm});
    hold on
    ylim([0 1]);
    xlim([0 40]);
    xlabel('eccentricity deg');
    ylabel('midget fraction');
    title('Our midget fraction');
    pbaspect([2 1 1]);
    drawnow
end
if p.Results.savePlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDir,'mRGCFractionByMeridian.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end

%% Plots related to midget RF density

figHandle = figure();
cardinalMeridianAngles=[0 90 180 270];
meridianColors={'r','b','g','k'};

for mm = 1:4
    
    meridianIdx = find(meridianAngles==cardinalMeridianAngles(mm),1);
    
    % load the empirical cone density measured by Curcio
    [coneDensitySqDeg, coneNativeSupportPosDeg] = getCurcioConeDensityByEccen(cardinalMeridianAngles(mm));
    % remove nan values
    isvalididx=find(~isnan(coneDensitySqDeg));
    coneNativeSupportPosDeg = coneNativeSupportPosDeg(isvalididx);
    coneDensitySqDeg = coneDensitySqDeg(isvalididx);
    
    % calculate the mRF density using Watson equation 8 at the sites of
    % empirical cone measurement
    [ mRFDensitySqDeg_watson ] = calcWatsonMidgetRFDensityByEccen(coneNativeSupportPosDeg, cardinalMeridianAngles(mm));
    
    % Plot the mRF density by eccentricity for Watson
    subplot(2,2,1);
    
    loglog(coneNativeSupportPosDeg(2:end),mRFDensitySqDeg_watson(2:end),'-','Color',meridianColors{mm});
    ylim([1e0 1e5]);
    xlabel('log10 eccentricity');
    ylabel('log10 mRF density / deg2');
    title('Watson''s mRF density by eccentricity');
    pbaspect([2 1 1]);
    hold on
    
    % Plot the mRF density by eccentricity from our functions
    subplot(2,2,2);
    mRFDensitySqDeg_ours = transformConeToMidgetRFDensity(coneDensityFit(coneNativeSupportPosDeg), ...
        'logitFitParams',fitParams(meridianIdx,1:2))';
    loglog(coneNativeSupportPosDeg(2:end),mRFDensitySqDeg_ours(2:end),'-','Color',meridianColors{mm});
    ylim([1e0 1e5]);
    xlabel('log10 eccentricity');
    ylabel('log10 mRF density / deg2');
    title('Our mRF density by eccentricity');
    pbaspect([2 1 1]);
    hold on
    
    % Plot the mRF : cone ratio for Watson
    subplot(2,2,3);
    loglog(coneNativeSupportPosDeg(2:end),mRFDensitySqDeg_watson(2:end)./coneDensitySqDeg(2:end),'-','Color',meridianColors{mm});
    ylim([1e-4 1e1]);
    xlabel('log10 eccentricity');
    ylabel('log10 mRF:cone');
    ylim([1e-4 1e2]);
    title('Watson''s mRF:cone ratio by eccentricity');
    pbaspect([2 1 1]);
    hold on
    
    % Plot the mRF : cone ratio for us
    subplot(2,2,4);
    mRFDensitySqDeg_ours = ...
        transformConeToMidgetRFDensity(coneDensityFit(coneNativeSupportPosDeg)','logitFitParams',fitParams(mm,1:2));
    loglog(coneNativeSupportPosDeg(2:end),mRFDensitySqDeg_ours(2:end)./coneDensitySqDeg(2:end),'-','Color',meridianColors{mm});
    ylim([1e-4 1e2]);
    xlabel('log10 eccentricity');
    ylabel('log10 mRF:cone');
    title('Our mRF:cone ratio by eccentricity');
    hold on
    pbaspect([2 1 1]);
    drawnow
    
end
if p.Results.savePlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDir,'mRFDensityAndRatioByMeridian.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end

end % function