function validation_MeridianPlots(varargin)
% Create plots that illustrate the behavior of the model along the cardinal
% merdians


%% Parse input and define variables
p = inputParser; p.KeepUnmatched = true;

% Optional anaysis params
p.addParameter('sampleResolutionDegreesRetina',0.01,@isnumeric);
p.addParameter('maxModeledEccentricityDegreesRetina',30,@isnumeric);
p.addParameter('meridianAngleResolutionDeg',90,@isnumeric);
p.addParameter('cardinalMeridianAngles',[0 90 180 270],@isnumeric);
p.addParameter('cardinalMeridianNames',{'nasal','superior','temporal','inferior'},@iscell);
p.addParameter('cardinalMeridianPlotColors',{'r','b','g','k'},@iscell);
p.addParameter('displacementMapPixelsPerDegRetina',10,@isnumeric);
p.addParameter('referenceEccenDegRetina',15,@isnumeric);

% Optional display and ouput params
p.addParameter('verbose',true,@islogical);
p.addParameter('savePlots',true,@islogical);
p.addParameter('pathToPlotOutputDirRoot','~/Desktop/rgcDisplacementMapPlots',@ischar);
p.addParameter('subjectName', '29986A', @ischar);

% parse
p.parse(varargin{:})


%% Setup
close all

% check or make a directory for output
if p.Results.savePlots
    if exist(fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName),'dir')==0
        mkdir(fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName));
    end
end

% Prepare the regular eccentricity support base
regularSupportPosDegRetina = ...
    0:p.Results.sampleResolutionDegreesRetina:p.Results.maxModeledEccentricityDegreesRetina;

% Define the cone and rgc density data to operate upon
coneDensityDataFileName = fullfile([getpref('rgcDisplacementMap','LocalDataPath') , '/Curcio_1990_JCompNeurol_HumanPhotoreceptorTopography/curcioRawConeDensity_',p.Results.subjectName,'.mat']);
rgcDensityDataFileName = fullfile([getpref('rgcDisplacementMap','LocalDataPath') , '/Curcio_1990_JCompNeurol_GanglionCellTopography/curcioRawRGCDensity_',p.Results.subjectName,'.mat']);

% make the displacement map for the cardinal meridians
[ ~, fitParams, ~, rgcDisplacementDegRetinaEachMeridian, mRGC_cumulativeEachMeridian, mRF_cumulativeEachMeridian ] = ...
    makeDisplacementMap(...
    'sampleResolutionDegreesRetina', p.Results.sampleResolutionDegreesRetina, ...
    'maxModeledEccentricityDegreesRetina', p.Results.maxModeledEccentricityDegreesRetina, ...
    'meridianAngleResolutionDeg', p.Results.meridianAngleResolutionDeg, ...
    'displacementMapPixelsPerDegRetina', p.Results.displacementMapPixelsPerDegRetina, ...
     'coneDensityDataFileName', coneDensityDataFileName, ...
     'rgcDensityDataFileName', rgcDensityDataFileName, ...
    'verbose', p.Results.verbose);

fitParams

% Plot the cumulatives and displacements across meridians
figHandle = figure();
for mm = 1:length(p.Results.cardinalMeridianAngles)
    
    % plot the displacement
    subplot(length(p.Results.cardinalMeridianAngles),2,mm*2);
    plot(regularSupportPosDegRetina,rgcDisplacementDegRetinaEachMeridian(mm,:),'-r')
    axis off;
    ylim([-.5 4.0]);
    if mm == length(p.Results.cardinalMeridianAngles)
        axis on;
        xlabel('eccentricity [deg retina]');
        ylabel('RGC displacement [deg retina]');
    end
    
    % Plot the cumulative functions
    subplot(length(p.Results.cardinalMeridianAngles),2,mm*2-1);
    plot(regularSupportPosDegRetina,mRGC_cumulativeEachMeridian(mm,:),'-k')
    axis off;
    if mm == length(p.Results.cardinalMeridianAngles)
        axis on;
        xlabel('eccentricity [deg retina]');
        ylabel('cells per sector');
    end
    hold on
    plot(regularSupportPosDegRetina,mRF_cumulativeEachMeridian(mm,:),'-b')
    ylim([0 8e5]);
    hold off
    drawnow
end
if p.Results.savePlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,'cumulativeAndDisplaceByMeridian.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end


% Show the cone --> mRF model
[ ~, figHandle ] = developMidgetRFFractionModel( 'makePlots', true );
hold on
% Add plot lines showing the fit by meridian
for mm = 1:length(p.Results.cardinalMeridianAngles)
    fitConeDensitySqDegRetina = getSplineFitToConeDensitySqDegRetina(p.Results.cardinalMeridianAngles(mm));
    coneDensitySqDegRetina = fitConeDensitySqDegRetina(regularSupportPosDegRetina);
    [ ~, mRFtoConeDensityRatio ] = transformConeToMidgetRFDensity( coneDensitySqDegRetina, 'linkingFuncParams', fitParams(mm,1:2) );
    xvals = log10(coneDensitySqDegRetina./max(coneDensitySqDegRetina));
    plot(xvals,mRFtoConeDensityRatio,'-','Color',p.Results.cardinalMeridianPlotColors{mm});
end
hold off
if p.Results.savePlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,'cone --> mRF RatioModel.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end

% Show the RGC --> mRGC model
[ ~, figHandle ] = developDaceyMidgetRGCFractionModel( 'makePlots', true );
hold on
% Add plot lines showing the fit by meridian
for mm = 1:length(p.Results.cardinalMeridianAngles)
    % Load the RGC Density Data from Curcio and Allen 1990
    [ RGCDensitySqDegRetina, nativeSupportPosDegRetina ] = loadRawRGCDensityByEccen( p.Results.cardinalMeridianAngles(mm) );
    
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
    plot(xvals,midgetFraction,'-','Color',p.Results.cardinalMeridianPlotColors{mm});
end
hold off
if p.Results.savePlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,'RGC --> mRGC RatioModel.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end

% Show the cone --> mRF model by merdian
figHandle = figure();
subplot(2,1,1)
for mm=1:4
    xFit= 1.4806e+04/100:100:1.4806e+04;
    tmp_mRFDensity = transformConeToMidgetRFDensity(xFit, 'linkingFuncParams', fitParams(mm,1:2));
    plot( log10(xFit ./ 1.4806e+04), tmp_mRFDensity./xFit,'-','Color',p.Results.cardinalMeridianPlotColors{mm});
    hold on
end
legend({p.Results.cardinalMeridianNames{:} 'fit'},'Location','southeast');
title('midget RF : cone ratio as a function of relative cone density');
xlabel('log10 proportion max cone density');
ylabel('mRF density : cone density');
legend(num2str(p.Results.cardinalMeridianAngles),'Location','southwest')
if p.Results.savePlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,'cone --> mRF model by meridian.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end



% Plot the spline cone density fits
figHandle = figure();
subplot(1,2,1)
for mm=1:4
    [coneDensitySqDegRetina, coneNativeSupportPosDegRetina] = loadRawConeDensityByEccen(p.Results.cardinalMeridianAngles(mm));
    dispConeNativeSupportPosDeg=coneNativeSupportPosDegRetina;
    dispConeNativeSupportPosDeg(1)=1e-2;
    [fitConeDensitySqDegRetina] = getSplineFitToConeDensitySqDegRetina(p.Results.cardinalMeridianAngles(mm));
    loglog(dispConeNativeSupportPosDeg,coneDensitySqDegRetina,'x','Color',p.Results.cardinalMeridianPlotColors{mm});
    xlim([1e-2,70]);
    hold on
    loglog(dispConeNativeSupportPosDeg,fitConeDensitySqDegRetina(coneNativeSupportPosDegRetina),'-','Color',p.Results.cardinalMeridianPlotColors{mm});
    xlabel('log10 Eccentricity [deg retina]');
    ylabel('log10 Cone density [counts / deg retina^2]');
end
legend(num2str(p.Results.cardinalMeridianAngles),'Location','southwest')

subplot(1,2,2)
interpolarMeridianAngles=[45 135 225 315];
for mm=1:4
    [fitConeDensitySqDegRetina] = getSplineFitToConeDensitySqDegRetina(interpolarMeridianAngles(mm));
    loglog(dispConeNativeSupportPosDeg,fitConeDensitySqDegRetina(coneNativeSupportPosDegRetina),'-','Color',p.Results.cardinalMeridianPlotColors{mm});
    xlim([1e-2,70]);
    xlabel('log10 Eccentricity [deg retina]');
    ylabel('log10 Cone density [counts / deg retina ^2]');
    hold on
end
legend(num2str(interpolarMeridianAngles),'Location','southwest')
if p.Results.savePlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,'coneDensitySplineModel.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end


% Plot the spline RGC density fits
figHandle = figure();
subplot(1,2,1)
for mm=1:4
    [RGCDensitySqDeg, RGCNativeSupportPosDeg] = loadRawRGCDensityByEccen(p.Results.cardinalMeridianAngles(mm));
    dispRGCNativeSupportPosDeg=RGCNativeSupportPosDeg;
    dispRGCNativeSupportPosDeg(1)=1e-2;
    [RGCDensityFit] = getSplineFitToRGCDensitySqDegRetina(p.Results.cardinalMeridianAngles(mm));
    loglog(dispRGCNativeSupportPosDeg,RGCDensitySqDeg,'x','Color',p.Results.cardinalMeridianPlotColors{mm});
    ylim([1,3000]);
    xlim([0.2,70]);
    hold on
    regularSupportPosDegRetina=1e-2:0.01:70;
    loglog(regularSupportPosDegRetina,RGCDensityFit(regularSupportPosDegRetina),'-','Color',p.Results.cardinalMeridianPlotColors{mm});
    xlabel('log10 Eccentricity [deg retina]');
    ylabel('log10 RGC density [counts / deg retina ^2]');
end
legend(num2str(p.Results.cardinalMeridianAngles),'Location','southwest')

subplot(1,2,2)
interpolarMeridianAngles=[45 135 225 315];
for mm=1:4
    [RGCDensityFit] = getSplineFitToRGCDensitySqDegRetina(interpolarMeridianAngles(mm));
    loglog(regularSupportPosDegRetina,RGCDensityFit(regularSupportPosDegRetina),'-','Color',p.Results.cardinalMeridianPlotColors{mm});
    xlim([0.2,70]);
    ylim([1,3000]);
    xlabel('log10 Eccentricity [deg retina]');
    ylabel('log10 RGC density [counts / deg retina ^2]');
    hold on
end
legend(num2str(interpolarMeridianAngles),'Location','southwest')
if p.Results.savePlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,'rgcDensitySplineModel.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end


%% Plot the spline cone and rgc fits for one meridian in linear space
mm=3;
figHandle = figure();
subplot(2,1,1)
[coneDensitySqDegRetina, coneNativeSupportPosDegRetina] = loadRawConeDensityByEccen(p.Results.cardinalMeridianAngles(mm));
dispConeNativeSupportPosDeg=coneNativeSupportPosDegRetina;
[fitConeDensitySqDegRetina] = getSplineFitToConeDensitySqDegRetina(p.Results.cardinalMeridianAngles(mm));
plot(dispConeNativeSupportPosDeg,coneDensitySqDegRetina,'x','Color',p.Results.cardinalMeridianPlotColors{mm});
xlim([0,30]);
ylim([0,3e4]);
hold on
regularSupportPosDegRetina=0:0.01:70;
plot(regularSupportPosDegRetina,fitConeDensitySqDegRetina(regularSupportPosDegRetina),'-','Color',p.Results.cardinalMeridianPlotColors{mm});
xlabel('Eccentricity [deg retina]');
ylabel('Cone density [counts / deg retina ^2]');
subplot(2,1,2)
[RGCDensitySqDeg, RGCNativeSupportPosDeg] = loadRawRGCDensityByEccen(p.Results.cardinalMeridianAngles(mm));
dispRGCNativeSupportPosDeg=RGCNativeSupportPosDeg;
[RGCDensityFit] = getSplineFitToRGCDensitySqDegRetina(p.Results.cardinalMeridianAngles(mm));
plot(dispRGCNativeSupportPosDeg,RGCDensitySqDeg,'x','Color',p.Results.cardinalMeridianPlotColors{mm});
xlim([0,30]);
ylim([0,2500]);
hold on
regularSupportPosDegRetina=1e-2:0.01:70;
loglog(regularSupportPosDegRetina,RGCDensityFit(regularSupportPosDegRetina),'-','Color',p.Results.cardinalMeridianPlotColors{mm});
xlabel('Eccentricity [deg]');
ylabel('RGC density [counts / deg retina ^2]');
if p.Results.savePlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,'splineFitConeRGCLinearSpace.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end

%% Plot the mRF and mRGC values for one meridian in linear space
mm=3;
figHandle = figure();
cardinalMeridianAngles=[0 90 180 270];
meridianColors={'r','b','g','k'};
subplot(2,1,1)
[fitConeDensitySqDegRetina] = getSplineFitToConeDensitySqDegRetina(p.Results.cardinalMeridianAngles(mm));
regularSupportPosDegRetina=0:0.01:70;
[ mRFDensitySqDeg ] = transformConeToMidgetRFDensity( fitConeDensitySqDegRetina(regularSupportPosDegRetina), 'linkingFuncParams', fitParams(mm,1:2) );
plot(regularSupportPosDegRetina,mRFDensitySqDeg,'-','Color',p.Results.cardinalMeridianPlotColors{mm});
xlim([0,30]);
ylim([0,3e4]);
xlabel('Eccentricity [deg retina]');
ylabel('mRF density [counts / deg retina ^2]');

subplot(2,1,2)
[RGCDensityFit] = getSplineFitToRGCDensitySqDegRetina(p.Results.cardinalMeridianAngles(mm));
[ mRGCDensitySqDeg ] = transformRGCToMidgetRGCDensityDacey( regularSupportPosDegRetina, RGCDensityFit(regularSupportPosDegRetina)', 'linkingFuncParams', fitParams(mm,3:end) );
plot(regularSupportPosDegRetina,mRGCDensitySqDeg,'-','Color',p.Results.cardinalMeridianPlotColors{mm});
xlim([0,30]);
ylim([0,2500]);
xlabel('Eccentricity [deg retina]');
ylabel('mRGC density [counts / deg retina^2]');
if p.Results.savePlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,'mRGCmRF_linearOneMeridian.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end


%% Plot the mRGC fraction for the cardinal meridians
figHandle = figure();
for mm = 1:4
        
    % Load the RGC Density Data from Curcio and Allen 1990:
    [ RGCDensitySqDeg, RGCNativeSupportPosDeg ] = loadRawRGCDensityByEccen( p.Results.cardinalMeridianAngles(mm) );
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
    midgetFraction_Dacey = calcDaceyMidgetFractionByEccenDegRetina(RGCNativeSupportPosDeg);
    plot(RGCNativeSupportPosDeg,midgetFraction_Dacey,'-k');
    xlabel('eccentricity deg (visual - Drasdo, retina - dacey)');
    ylabel('midget fraction');
    ylim([0 1]);
    xlim([0 40]);
    title('Dacey and Drasdo midget fraction');
    pbaspect([2 1 1]);
    
    % Plot our midget fraction
    subplot(1,2,2);
    [ ~, midgetFraction_ours ] = transformRGCToMidgetRGCDensityDacey( RGCNativeSupportPosDeg, RGCDensitySqDeg, 'linkingFuncParams', fitParams(mm,3:end) );
    plot(RGCNativeSupportPosDeg,midgetFraction_ours,'-','Color',p.Results.cardinalMeridianPlotColors{mm});
    hold on
    ylim([0 1]);
    xlim([0 40]);
    xlabel('eccentricity deg retina');
    ylabel('midget fraction');
    title('Our midget fraction');
    pbaspect([2 1 1]);
    drawnow
end
if p.Results.savePlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,'mRGCFractionByMeridian.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end

%% Plots related to midget RF density

figHandle = figure();
for mm = 1:4
        
    % Load the Curcio average data reported in the paper
    [coneDensitySqDegRetina, coneDensitySupportDegRetina] = loadRawConeDensityByEccen(p.Results.cardinalMeridianAngles(mm));

    % Convert to units of visual degrees
    coneDensitySupportMmRetina = convert_degRetina_to_mmRetina(coneDensitySupportDegRetina);
    coneDensitySupportDegVisual = convert_mmRetina_to_degVisual(coneDensitySupportMmRetina, p.Results.cardinalMeridianAngles(mm));
    coneDensitySqMmRetina = coneDensitySqDegRetina .* calc_degSqRetina_per_mmSqRetina();
    coneDensitySqDegVisual = coneDensitySqMmRetina ./ calc_degSqVisual_per_mmSqRetina(coneDensitySupportMmRetina, p.Results.cardinalMeridianAngles(mm));
    
    % calculate the mRF density using Watson equation 8 at the sites of
    % empirical cone measurement
    [ mRFDensitySqDegVisual_watson ] = calcWatsonMidgetRFDensityByEccenDegVisual(coneDensitySupportDegVisual, p.Results.cardinalMeridianAngles(mm));
    
    % Plot the mRF density by eccentricity for Watson
    subplot(2,2,1);
    
    loglog(coneDensitySupportDegVisual(2:end),mRFDensitySqDegVisual_watson(2:end),'-','Color',p.Results.cardinalMeridianPlotColors{mm});
    ylim([1e0 1e5]);
    xlim([1e-1 1e2]);
    xlabel('log10 eccentricity deg visual');
    ylabel('log10 mRF density / deg visual ^2');
    title('Watson''s mRF density by eccentricity');
    pbaspect([2 1 1]);
    hold on
    
    % Plot the mRF density by eccentricity from our functions
    subplot(2,2,2);
    mRFDensitySqDegRetina_ours = transformConeToMidgetRFDensity(coneDensitySqDegRetina, ...
        'linkingFuncParams',fitParams(mm,1:2));
    % convert from SqDegRetina to SqDegVisual field via SqMmRetina
    mRFDensitySqMmRetina_ours = mRFDensitySqDegRetina_ours .* calc_degSqRetina_per_mmSqRetina();
    mRFDensitySqDegVisual_ours = mRFDensitySqMmRetina_ours ./ calc_degSqVisual_per_mmSqRetina(coneDensitySupportMmRetina, p.Results.cardinalMeridianAngles(mm));
    
    loglog(coneDensitySupportDegVisual(2:end),mRFDensitySqDegVisual_ours(2:end),'-','Color',p.Results.cardinalMeridianPlotColors{mm});
    ylim([1e0 1e5]);
    xlim([1e-1 1e2]);
    xlabel('log10 eccentricity visual deg');
    ylabel('log10 mRF density / visual deg2');
    title('Our mRF density by visual eccentricity');
    pbaspect([2 1 1]);
    hold on
    
    % Plot the mRF : cone ratio for Watson
    subplot(2,2,3);
    loglog(coneDensitySupportDegVisual(2:end),mRFDensitySqDegVisual_watson(2:end)./coneDensitySqDegVisual(2:end),'-','Color',p.Results.cardinalMeridianPlotColors{mm});
    ylim([1e-2 10]);
    xlim([1e-1 1e2]);
    xlabel('log10 eccentricity visual degrees');
    ylabel('log10 mRF:cone');
    title('Watson''s mRF:cone ratio by visual eccentricity');
    pbaspect([2 1 1]);
    hold on
    
    
    % Plot the mRF : cone ratio for us
    subplot(2,2,4);
    loglog(coneDensitySupportDegVisual(2:end),mRFDensitySqDegVisual_ours(2:end)./coneDensitySqDegVisual(2:end),'-','Color',p.Results.cardinalMeridianPlotColors{mm});
    ylim([1e-2 10]);
    xlim([1e-1 1e2]);
    xlabel('log10 eccentricity visual degrees');
    ylabel('log10 mRF:cone');
    title('Our mRF:cone ratio by visual eccentricity');
    hold on
    pbaspect([2 1 1]);
    drawnow
    
end
if p.Results.savePlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,'mRFDensityAndRatioByMeridian.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end

end % function