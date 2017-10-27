function validation_MeridianPlots(varargin)
% Create plots that illustrate the behavior of the model along the cardinal
% merdians


%% Parse input and define variables
p = inputParser; p.KeepUnmatched = true;

% Optional anaysis params
p.addParameter('sampleResolutionDegrees',0.01,@isnumeric);
p.addParameter('maxModeledEccentricity',30,@isnumeric);
p.addParameter('meridianAngleResolutionDeg',90,@isnumeric);
p.addParameter('meridianNames',{'nasal','superior','temporal','inferior'},@iscell);
p.addParameter('displacementMapPixelsPerDeg',10,@isnumeric);
p.addParameter('referenceEccenDegRetina',15,@isnumeric);

% Optional display and ouput params
p.addParameter('verbose',true,@islogical);
p.addParameter('savePlots',true,@islogical);
p.addParameter('pathToPlotOutputDir','~/Desktop/rgcDisplacementMapPlots',@ischar);

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
regularSupportPosDegRetina = ...
    0:p.Results.sampleResolutionDegrees:p.Results.maxModeledEccentricity;

% make the displacement map for the cardinal meridians
[ ~, fitParams, meridianAngles, rgcDisplacementEachMeridian, mRGC_cumulativeEachMeridian, mRF_cumulativeEachMeridian ] = ...
    makeDisplacementMap(...
    'sampleResolutionDegrees', p.Results.sampleResolutionDegrees, ...
    'maxModeledEccentricity', p.Results.maxModeledEccentricity, ...
    'meridianAngleResolutionDeg', p.Results.meridianAngleResolutionDeg, ...
    'displacementMapPixelsPerDeg', p.Results.displacementMapPixelsPerDeg, ...
    'verbose', p.Results.verbose);

fitParams

% Plot the cumulatives and displacements across meridians
figHandle = figure();
for mm = 1:length(meridianAngles)
    
    % plot the displacement
    subplot(length(meridianAngles),2,mm*2);
    plot(regularSupportPosDegRetina,rgcDisplacementEachMeridian(mm,:),'-r')
    axis off;
    ylim([-.5 4.0]);
    if mm == length(meridianAngles)
        axis on;
        xlabel('eccentricity [deg]');
        ylabel('RGC displacement [deg]');
    end
    
    % Plot the cumulative functions
    subplot(length(meridianAngles),2,mm*2-1);
    plot(regularSupportPosDegRetina,mRGC_cumulativeEachMeridian(mm,:),'-k')
    axis off;
    if mm == length(meridianAngles)
        axis on;
        xlabel('eccentricity [deg]');
        ylabel('cells per sector');
    end
    hold on
    plot(regularSupportPosDegRetina,mRF_cumulativeEachMeridian(mm,:),'-b')
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
hold on
% Add plot lines showing the fit by meridian
meridianColors={'r','b','g','k'};
for mm = 1:length(meridianAngles)
    coneDensityFit = getSplineFitToConeDensity(meridianAngles(mm));
    coneDensitySqDeg = coneDensityFit(regularSupportPosDegRetina);
    [ ~, mRFtoConeDensityRatio ] = transformConeToMidgetRFDensity( coneDensitySqDeg, 'linkingFuncParams', fitParams(mm,1:2) );
    xvals = log10(coneDensitySqDeg./max(coneDensitySqDeg));
    plot(xvals,mRFtoConeDensityRatio,'-','Color',meridianColors{mm});
end
hold off
if p.Results.savePlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDir,'cone --> mRF RatioModel.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end

% Show the RGC --> mRGC model
[ ~, figHandle ] = developDaceyMidgetRGCFractionModel( 'makePlots', true );
hold on
% Add plot lines showing the fit by meridian
meridianColors={'r','b','g','k'};
for mm = 1:length(meridianAngles)
    % Load the RGC Density Data from Curcio and Allen 1990
    [ RGCDensitySqDegRetina, nativeSupportPosDegRetina ] = loadRawRGCDensityByEccen( meridianAngles(mm) );
    
    % remove nan values
    isvalididx=find(~isnan(RGCDensitySqDegRetina)  );
    nativeSupportPosDegRetina = nativeSupportPosDegRetina(isvalididx);
    RGCDensitySqDegRetina = RGCDensitySqDegRetina(isvalididx);
    
    % Fit a spline to the RGC density data
    RGCDensitySqDegRetinaFit = fit(nativeSupportPosDegRetina',RGCDensitySqDegRetina','smoothingspline', 'Exclude',find(isnan(RGCDensitySqDegRetina)),'SmoothingParam', 1);
    
    % Obtain the cumulative RGC function
    RGC_ringcount = calcCumulative(regularSupportPosDegRetina,RGCDensitySqDegRetinaFit(regularSupportPosDegRetina)');
    
    % Find the index position in the regularSupportPosDegRetina that is as close
    % as possible to the referenceEccenDegRetinaDegRetina
    refPointIdx= ...
        find((regularSupportPosDegRetina-p.Results.referenceEccenDegRetina)== ...
        min(abs(regularSupportPosDegRetina-p.Results.referenceEccenDegRetina)));

    % Calculate a proportion of the cumulative RGC density counts, relative
    % to the reference point (which is assigned a value of unity)
    propRGC_ringcount=RGC_ringcount./RGC_ringcount(refPointIdx);

    [ ~, midgetFraction ] = transformRGCToMidgetRGCDensityDacey( regularSupportPosDegRetina, RGCDensitySqDegRetinaFit(regularSupportPosDegRetina)', 'linkingFuncParams', fitParams(mm,3:end) );
    xvals = propRGC_ringcount;
    plot(xvals,midgetFraction,'-','Color',meridianColors{mm});
end
hold off
if p.Results.savePlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDir,'RGC --> mRGC RatioModel.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end

% Show the cone --> mRF model by merdian
figHandle = figure();
cardinalMeridianAngles=[0 90 180 270];
meridianColors={'r','b','g','k'};
subplot(2,1,1)
for mm=1:4
    xFit= 1.4806e+04/100:100:1.4806e+04;
    tmp_mRFDensity = transformConeToMidgetRFDensity(xFit, 'linkingFuncParams', fitParams(mm,1:2));
    plot( log10(xFit ./ 1.4806e+04), tmp_mRFDensity./xFit,'-','Color',meridianColors{mm});
    hold on
end
legend({p.Results.meridianNames{:} 'fit'},'Location','southeast');
title('midget RF : cone ratio as a function of relative cone density');
xlabel('log10 proportion max cone density');
ylabel('mRF density : cone density');
legend(num2str(cardinalMeridianAngles),'Location','southwest')
if p.Results.savePlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDir,'cone --> mRF model by meridian.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end



% Plot the spline cone density fits
figHandle = figure();
cardinalMeridianAngles=[0 90 180 270];
meridianColors={'r','b','g','k'};
subplot(1,2,1)
for mm=1:4
    [coneDensitySqDeg, coneNativeSupportPosDeg] = loadRawConeDensityByEccen(cardinalMeridianAngles(mm));
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
    [RGCDensitySqDeg, RGCNativeSupportPosDeg] = loadRawRGCDensityByEccen(cardinalMeridianAngles(mm));
    dispRGCNativeSupportPosDeg=RGCNativeSupportPosDeg;
    dispRGCNativeSupportPosDeg(1)=1e-2;
    [RGCDensityFit] = getSplineFitToRGCDensity(cardinalMeridianAngles(mm));
    loglog(dispRGCNativeSupportPosDeg,RGCDensitySqDeg,'x','Color',meridianColors{mm});
    xlim([1e-2,70]);
    hold on
    regularSupportPosDegRetina=1e-2:0.01:70;
    loglog(regularSupportPosDegRetina,RGCDensityFit(regularSupportPosDegRetina),'-','Color',meridianColors{mm});
    xlabel('log10 Eccentricity [deg]');
    ylabel('log10 RGC density [counts / deg2]');
end
legend(num2str(cardinalMeridianAngles),'Location','southwest')

subplot(1,2,2)
interpolarMeridianAngles=[45 135 225 315];
for mm=1:4
    [RGCDensityFit] = getSplineFitToRGCDensity(interpolarMeridianAngles(mm));
    loglog(regularSupportPosDegRetina,RGCDensityFit(regularSupportPosDegRetina),'-','Color',meridianColors{mm});
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


%% Plot the spline cone and rgc fits for one meridian in linear space
mm=3;
figHandle = figure();
cardinalMeridianAngles=[0 90 180 270];
meridianColors={'r','b','g','k'};
subplot(2,1,1)
[coneDensitySqDeg, coneNativeSupportPosDeg] = loadRawConeDensityByEccen(cardinalMeridianAngles(mm));
dispConeNativeSupportPosDeg=coneNativeSupportPosDeg;
[coneDensityFit] = getSplineFitToConeDensity(cardinalMeridianAngles(mm));
plot(dispConeNativeSupportPosDeg,coneDensitySqDeg,'x','Color',meridianColors{mm});
xlim([0,30]);
ylim([0,3e4]);
hold on
regularSupportPosDegRetina=0:0.01:70;
plot(regularSupportPosDegRetina,coneDensityFit(regularSupportPosDegRetina),'-','Color',meridianColors{mm});
xlabel('Eccentricity [deg]');
ylabel('Cone density [counts / deg2]');
subplot(2,1,2)
[RGCDensitySqDeg, RGCNativeSupportPosDeg] = loadRawRGCDensityByEccen(cardinalMeridianAngles(mm));
dispRGCNativeSupportPosDeg=RGCNativeSupportPosDeg;
[RGCDensityFit] = getSplineFitToRGCDensity(cardinalMeridianAngles(mm));
plot(dispRGCNativeSupportPosDeg,RGCDensitySqDeg,'x','Color',meridianColors{mm});
xlim([0,30]);
ylim([0,2500]);
hold on
regularSupportPosDegRetina=1e-2:0.01:70;
loglog(regularSupportPosDegRetina,RGCDensityFit(regularSupportPosDegRetina),'-','Color',meridianColors{mm});
xlabel('Eccentricity [deg]');
ylabel('RGC density [counts / deg2]');
if p.Results.savePlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDir,'splineFitConeRGCLinearSpace.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end

%% Plot the mRF and mRGC values for one meridian in linear space
mm=3;
figHandle = figure();
cardinalMeridianAngles=[0 90 180 270];
meridianColors={'r','b','g','k'};
subplot(2,1,1)
[coneDensityFit] = getSplineFitToConeDensity(cardinalMeridianAngles(mm));
regularSupportPosDegRetina=0:0.01:70;
[ mRFDensitySqDeg ] = transformConeToMidgetRFDensity( coneDensityFit(regularSupportPosDegRetina), 'linkingFuncParams', fitParams(mm,1:2) );
plot(regularSupportPosDegRetina,mRFDensitySqDeg,'-','Color',meridianColors{mm});
xlim([0,30]);
ylim([0,3e4]);
xlabel('Eccentricity [deg]');
ylabel('mRF density [counts / deg2]');

subplot(2,1,2)
[RGCDensityFit] = getSplineFitToRGCDensity(cardinalMeridianAngles(mm));
[ mRGCDensitySqDeg ] = transformRGCToMidgetRGCDensityDacey( regularSupportPosDegRetina, RGCDensityFit(regularSupportPosDegRetina)', 'linkingFuncParams', fitParams(mm,3:end) );
plot(regularSupportPosDegRetina,mRGCDensitySqDeg,'-','Color',meridianColors{mm});
xlim([0,30]);
ylim([0,2500]);
xlabel('Eccentricity [deg]');
ylabel('mRGC density [counts / deg2]');
if p.Results.savePlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDir,'mRGCmRF_linearOneMeridian.pdf');
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
    [ RGCDensitySqDeg, RGCNativeSupportPosDeg ] = loadRawRGCDensityByEccen( cardinalMeridianAngles(mm) );
    % remove nan values
    isvalididx=find(~isnan(RGCDensitySqDeg)  );
    RGCNativeSupportPosDeg = RGCNativeSupportPosDeg(isvalididx);
    RGCDensitySqDeg = RGCDensitySqDeg(isvalididx);
    
    % Plot the Drasdo and Dacey midget fraction
    subplot(1,2,1);
    f0 = 0.8928; rm = 41.03; % Drasdo's values
    midgetFraction_Drasdo = calcDrasdoMidgetFractionByVisualEccen(RGCNativeSupportPosDeg,f0,rm);
    plot(RGCNativeSupportPosDeg,midgetFraction_Drasdo,'-k');
    hold on
    midgetFraction_Dacey = calcDaceyMidgetFractionByVisualEccen(RGCNativeSupportPosDeg);
    plot(RGCNativeSupportPosDeg,midgetFraction_Dacey,'-k');
    xlabel('eccentricity deg (visual - Drasdo, retina - dacey)');
    ylabel('midget fraction');
    ylim([0 1]);
    xlim([0 40]);
    title('Dacey and Drasdo midget fraction');
    pbaspect([2 1 1]);
    
    % Plot our midget fraction
    subplot(1,2,2);
    [ ~, midgetFraction_ours ] = transformRGCToMidgetRGCDensityDacey( RGCNativeSupportPosDeg, RGCDensitySqDeg, 'linkingFuncParams', fitParams(meridianIdx,3:end) );
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
    [coneDensitySqDeg, coneNativeSupportPosDeg] = loadRawConeDensityByEccen(cardinalMeridianAngles(mm));
    % remove nan values
    isvalididx=find(~isnan(coneDensitySqDeg));
    coneNativeSupportPosDeg = coneNativeSupportPosDeg(isvalididx);
    coneDensitySqDeg = coneDensitySqDeg(isvalididx);

    % Get the spline fit to cone density for this medidain
    [coneDensityFit] = getSplineFitToConeDensity(cardinalMeridianAngles(mm));

    % calculate the mRF density using Watson equation 8 at the sites of
    % empirical cone measurement
    [ mRFDensitySqDeg_watson ] = calcWatsonMidgetRFDensityByEccen(coneNativeSupportPosDeg, cardinalMeridianAngles(mm));
    
    % Plot the mRF density by eccentricity for Watson
    subplot(2,2,1);
    
    loglog(coneNativeSupportPosDeg(2:end),mRFDensitySqDeg_watson(2:end),'-','Color',meridianColors{mm});
    ylim([1e0 1e5]);
    xlim([1e-1 1e2]);
    xlabel('log10 eccentricity');
    ylabel('log10 mRF density / deg2');
    title('Watson''s mRF density by eccentricity');
    pbaspect([2 1 1]);
    hold on
    
    % Plot the mRF density by eccentricity from our functions
    subplot(2,2,2);
    mRFDensitySqDeg_ours = transformConeToMidgetRFDensity(coneDensityFit(coneNativeSupportPosDeg), ...
        'linkingFuncParams',fitParams(meridianIdx,1:2))';
    loglog(coneNativeSupportPosDeg(2:end),mRFDensitySqDeg_ours(2:end),'-','Color',meridianColors{mm});
    ylim([1e0 1e5]);
    xlim([1e-1 1e2]);
    xlabel('log10 eccentricity');
    ylabel('log10 mRF density / deg2');
    title('Our mRF density by eccentricity');
    pbaspect([2 1 1]);
    hold on
    
    % Plot the mRF : cone ratio for Watson
    subplot(2,2,3);
    loglog(coneNativeSupportPosDeg(2:end),mRFDensitySqDeg_watson(2:end)./coneDensitySqDeg(2:end),'-','Color',meridianColors{mm});
    ylim([1e-2 10]);
    xlim([1e-1 1e2]);
    xlabel('log10 eccentricity');
    ylabel('log10 mRF:cone');
    title('Watson''s mRF:cone ratio by eccentricity');
    pbaspect([2 1 1]);
    hold on
    
    % Plot the mRF : cone ratio for us
    subplot(2,2,4);
    mRFDensitySqDeg_ours = ...
        transformConeToMidgetRFDensity(coneDensityFit(coneNativeSupportPosDeg)','linkingFuncParams',fitParams(mm,1:2));
    loglog(coneNativeSupportPosDeg(2:end),mRFDensitySqDeg_ours(2:end)./coneDensitySqDeg(2:end),'-','Color',meridianColors{mm});
    ylim([1e-2 10]);
    xlim([1e-1 1e2]);
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