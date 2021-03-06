function demo_MeridianPlots(varargin)
% Demonstrate model performance on the cardinal meridians
%
% Description:
%   This validation function runs the displacement model for each of the
%   cardinal meridians and then saves a set of plots that illustrate the
%   results. The analysis parameters are explicitly defined here and then
%   passed to the main, createDisplacementModel routine.
%
%   The "subjectName" corresponds to one of the Curcio 1990 datasets
%   present within the data directory of this toolbox. Typical options
%   include:
%       reportedAverage   - The values reported in the Curcio 1990 papers
%       computedAverage   - Our derivation of average values from the 
%                           Curcio data
%       29986A            - Data from subject 29986, averaged over both
%                           eyes
%


%% Parse input and define variables
p = inputParser; p.KeepUnmatched = true;

% Optional anaysis params
p.addParameter('sampleResolutionDegVisual',0.01,@isnumeric);
p.addParameter('maxModeledEccentricityDegVisual',40,@isnumeric);
p.addParameter('meridianAngleResolutionDeg',90,@isnumeric);
p.addParameter('cardinalMeridianAngles',[0 90 180 270],@isnumeric);
p.addParameter('cardinalMeridianNames',{'nasal','superior','temporal','inferior'},@iscell);
p.addParameter('cardinalMeridianPlotColors',{'r','b','g','k'},@iscell);
p.addParameter('displacementMapPixelsPerDegVisual',10,@isnumeric);
p.addParameter('referenceEccenDegVisual',10,@isnumeric);
p.addParameter('subjectName', 'computedAverage', @ischar);

% Optional display and ouput params
p.addParameter('verbose',true,@islogical);
p.addParameter('saveVectorPlots',true,@islogical);
p.addParameter('pathToPlotOutputDirRoot','~/Desktop/rgcDisplacementMapPlots',@ischar);

% parse
p.parse(varargin{:})


%% Setup
close all

% check or make a directory for output
if p.Results.saveVectorPlots
    if exist(fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName),'dir')==0
        mkdir(fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName));
    end
end

% Define the cone and rgc density data to operate upon
coneDensityDataFileName = fullfile([getpref('rgcDisplacementMap','LocalDataPath') , '/Curcio_1990_JCompNeurol_HumanPhotoreceptorTopography/curcioRawConeDensity_',p.Results.subjectName,'.mat']);
rgcDensityDataFileName = fullfile([getpref('rgcDisplacementMap','LocalDataPath') , '/Curcio_1990_JCompNeurol_GanglionCellTopography/curcioRawRGCDensity_',p.Results.subjectName,'.mat']);

% Create the displacement model
[ rgcDisplacementByMeridian, meridianAngleSupport, regularSupportPosDegVisual, ~, mRF_RingCumulativeByMeridian, mRGC_RingCumulativeByMeridian, fitParamsByMeridian, ~, ~ ] = ...
    createDisplacementModel(...
    'sampleResolutionDegVisual', p.Results.sampleResolutionDegVisual, ...
    'maxModeledEccentricityDegVisual', p.Results.maxModeledEccentricityDegVisual, ...
    'meridianAngleResolutionDeg', p.Results.meridianAngleResolutionDeg, ...
    'displacementMapPixelsPerDegVisual', p.Results.displacementMapPixelsPerDegVisual, ...
    'coneDensityDataFileName', coneDensityDataFileName, ...
    'rgcDensityDataFileName', rgcDensityDataFileName, ...
    'verbose', p.Results.verbose);

% Plot the cumulatives and displacements across meridians
figHandle = figure();
for mm = 1:length(p.Results.cardinalMeridianAngles)
    
    % plot the displacement
    subplot(length(p.Results.cardinalMeridianAngles),2,mm*2);
    plot(regularSupportPosDegVisual,rgcDisplacementByMeridian(mm,:),'-r')
    axis off;
    ylim([-.5 4.0]);
    title(p.Results.cardinalMeridianNames{mm});
    if mm == length(p.Results.cardinalMeridianAngles)
        axis on;
        xlabel('eccentricity [deg retina]');
        ylabel('RGC displacement [deg retina]');
    end
    
    % Plot the cumulative functions
    subplot(length(p.Results.cardinalMeridianAngles),2,mm*2-1);
    plot(regularSupportPosDegVisual,mRGC_RingCumulativeByMeridian(mm,:),'-k')
    axis off;
    title(p.Results.cardinalMeridianNames{mm});
    if mm == length(p.Results.cardinalMeridianAngles)
        axis on;
        xlabel('eccentricity [deg retina]');
        ylabel('cells per sector');
    end
    hold on
    plot(regularSupportPosDegVisual,mRF_RingCumulativeByMeridian(mm,:),'-b')
    ylim([0 8e5]);
    hold off
    drawnow
    
end
if p.Results.saveVectorPlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,'cumulativeAndDisplaceByMeridian.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end


% Test the displacement
figHandle = figure();
hold on
for mm = 1:length(p.Results.cardinalMeridianAngles)
    % Load the cone density function
    fitConeDensitySqDegVisual = getSplineFitToConeDensitySqDegVisual(p.Results.cardinalMeridianAngles(mm));
    coneDensityOverRegularSupport = ...
        zeroOpticDiscPoints(fitConeDensitySqDegVisual(regularSupportPosDegVisual),regularSupportPosDegVisual, meridianAngleSupport(mm));
    [ mRFDensitySqDegVisual ] = transformConeToMidgetRFDensity( coneDensityOverRegularSupport, 'linkingFuncParams', fitParamsByMeridian(mm,1:2) );
    
    % Load the RGC density function
    fitRGCDensitySqDegVisual = getSplineFitToRGCDensitySqDegVisual(p.Results.cardinalMeridianAngles(mm));
    RGCDensityOverRegularSupport = ...
        zeroOpticDiscPoints(fitRGCDensitySqDegVisual(regularSupportPosDegVisual),regularSupportPosDegVisual, meridianAngleSupport(mm));
    [ mRGCDensitySqDegVisual ] = transformRGCToMidgetRGCDensityDacey( regularSupportPosDegVisual, RGCDensityOverRegularSupport, 'linkingFuncParams', fitParamsByMeridian(mm,3:end) );
    extendedRegularSupportPosDegVisual = [regularSupportPosDegVisual (regularSupportPosDegVisual(end)+diff(regularSupportPosDegVisual(1:2)))];
    ringArea = diff(extendedRegularSupportPosDegVisual.^2 * pi);
    mRGCCountsPerSector = mRGCDensitySqDegVisual .* ringArea;
    mRFCountsPerSector = mRFDensitySqDegVisual .* ringArea;
    
    % Displace the mRGC density values to the cone locations
    mRGCCountsPerSectorAtConeLocations = applyDisplacement(mRGCCountsPerSector, rgcDisplacementByMeridian(mm,:), regularSupportPosDegVisual);
    
    % Down-sample the displaced values and support
    binSize = 25;
    coarse_regularSupportPosDegVisual=[];
    coarse_mRGCCountsPerSectorAtConeLocations=[];
    for kk = 1:floor(length(regularSupportPosDegVisual)/binSize)
        coarse_regularSupportPosDegVisual(kk)=mean(regularSupportPosDegVisual((kk-1)*binSize+1:kk*binSize));
        coarse_mRGCCountsPerSectorAtConeLocations(kk)=mean(mRGCCountsPerSectorAtConeLocations((kk-1)*binSize+1:kk*binSize));
    end
    
    % plot the displacement
    subplot(length(p.Results.cardinalMeridianAngles),1,mm);
    plot(regularSupportPosDegVisual,mRFCountsPerSector,'-k')
    hold on
    plot(coarse_regularSupportPosDegVisual,coarse_mRGCCountsPerSectorAtConeLocations,'.r')
    title(p.Results.cardinalMeridianNames{mm});
    if mm == length(p.Results.cardinalMeridianAngles)
        axis on;
        xlabel('eccentricity [deg retina]');
        ylabel('counts per sector');
        legend({'mRF','displaced mRGC'});
    end
end
if p.Results.saveVectorPlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,'mRGCDensityDisplacedToConeLocations.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end



% Show the cone --> mRF model
[ ~, figHandle ] = developMidgetRFFractionModel( 'makePlots', true );
hold on
% Add plot lines showing the fit by meridian
for mm = 1:length(p.Results.cardinalMeridianAngles)
    fitConeDensitySqDegVisual = getSplineFitToConeDensitySqDegVisual(p.Results.cardinalMeridianAngles(mm));
    coneDensityOverRegularSupport = ...
        zeroOpticDiscPoints(fitConeDensitySqDegVisual(regularSupportPosDegVisual),regularSupportPosDegVisual, meridianAngleSupport(mm));
    [ ~, mRFtoConeDensityRatio ] = transformConeToMidgetRFDensity( coneDensityOverRegularSupport, 'linkingFuncParams', fitParamsByMeridian(mm,1:2) );
    xvals = log10(coneDensityOverRegularSupport./max(coneDensityOverRegularSupport));
    plot(xvals,mRFtoConeDensityRatio,'-','Color',p.Results.cardinalMeridianPlotColors{mm});
end
hold off
if p.Results.saveVectorPlots
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
    [ RGCDensitySqDegVisual, nativeSupportPosDegVisual ] = loadRawRGCDensityByEccen( p.Results.cardinalMeridianAngles(mm) );
    
    % remove nan values
    isvalididx=find(~isnan(RGCDensitySqDegVisual)  );
    nativeSupportPosDegVisual = nativeSupportPosDegVisual(isvalididx);
    RGCDensitySqDegVisual = RGCDensitySqDegVisual(isvalididx);
    
    % Fit a spline to the RGC density data
    RGCDensitySqDegVisualFit = fit(nativeSupportPosDegVisual',RGCDensitySqDegVisual','smoothingspline', 'Exclude',find(isnan(RGCDensitySqDegVisual)),'SmoothingParam', 1);
    RGCDensityOverRegularSupport = ...
        zeroOpticDiscPoints(RGCDensitySqDegVisualFit(regularSupportPosDegVisual),regularSupportPosDegVisual, meridianAngleSupport(mm));

    % Obtain the cumulative RGC function
    RGC_ringcount = calcRingCumulative(regularSupportPosDegVisual,RGCDensityOverRegularSupport');
    
    % Find the index position in the regularSupportPosDegVisual that is as close
    % as possible to the referenceEccenDegVisualDegVisual
    refPointIdx= ...
        find((regularSupportPosDegVisual-p.Results.referenceEccenDegVisual)== ...
        min(abs(regularSupportPosDegVisual-p.Results.referenceEccenDegVisual)));
    
    % Calculate a proportion of the cumulative RGC density counts, relative
    % to the reference point (which is assigned a value of unity)
    propRGC_ringcount=RGC_ringcount./RGC_ringcount(refPointIdx);
    
    [ ~, midgetFraction ] = transformRGCToMidgetRGCDensityDacey( regularSupportPosDegVisual, RGCDensityOverRegularSupport', 'linkingFuncParams', fitParamsByMeridian(mm,3:end) );
    xvals = propRGC_ringcount;
    plot(xvals,midgetFraction,'-','Color',p.Results.cardinalMeridianPlotColors{mm});
end
hold off
if p.Results.saveVectorPlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,'RGC --> mRGC RatioModel.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end

% Show the cone --> mRF model by merdian
figHandle = figure();
subplot(2,1,1)
for mm=1:4
    xFit= 1.4806e+04/100:100:1.4806e+04;
    tmp_mRFDensity = transformConeToMidgetRFDensity(xFit, 'linkingFuncParams', fitParamsByMeridian(mm,1:2));
    plot( log10(xFit ./ 1.4806e+04), tmp_mRFDensity./xFit,'-','Color',p.Results.cardinalMeridianPlotColors{mm});
    hold on
end
legend({p.Results.cardinalMeridianNames{:} 'fit'},'Location','southeast');
title('midget RF : cone ratio as a function of relative cone density');
xlabel('log10 proportion max cone density');
ylabel('mRF density : cone density');
legend(num2str(p.Results.cardinalMeridianAngles),'Location','southwest')
if p.Results.saveVectorPlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,'cone --> mRF model by meridian.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end



% Plot the spline cone density fits
figHandle = figure();
subplot(1,2,1)
for mm=1:4
    [coneDensitySqDegVisual, coneNativeSupportPosDegVisual] = loadRawConeDensityByEccen(p.Results.cardinalMeridianAngles(mm));
    dispConeNativeSupportPosDeg=coneNativeSupportPosDegVisual;
    dispConeNativeSupportPosDeg(1)=1e-2;
    [fitConeDensitySqDegVisual] = getSplineFitToConeDensitySqDegVisual(p.Results.cardinalMeridianAngles(mm));
    loglog(dispConeNativeSupportPosDeg,coneDensitySqDegVisual,'x','Color',p.Results.cardinalMeridianPlotColors{mm});
    xlim([1e-2,70]);
    hold on
    loglog(dispConeNativeSupportPosDeg,fitConeDensitySqDegVisual(coneNativeSupportPosDegVisual),'-','Color',p.Results.cardinalMeridianPlotColors{mm});
    xlabel('log10 Eccentricity [deg retina]');
    ylabel('log10 Cone density [counts / deg retina^2]');
end
legend(num2str(p.Results.cardinalMeridianAngles),'Location','southwest')

subplot(1,2,2)
interpolarMeridianAngles=[45 135 225 315];
for mm=1:4
    [fitConeDensitySqDegVisual] = getSplineFitToConeDensitySqDegVisual(interpolarMeridianAngles(mm));
    loglog(dispConeNativeSupportPosDeg,fitConeDensitySqDegVisual(coneNativeSupportPosDegVisual),'-','Color',p.Results.cardinalMeridianPlotColors{mm});
    xlim([1e-2,70]);
    xlabel('log10 Eccentricity [deg retina]');
    ylabel('log10 Cone density [counts / deg retina ^2]');
    hold on
end
legend(num2str(interpolarMeridianAngles),'Location','southwest')
if p.Results.saveVectorPlots
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
    [RGCDensityFit] = getSplineFitToRGCDensitySqDegVisual(p.Results.cardinalMeridianAngles(mm));
    loglog(dispRGCNativeSupportPosDeg,RGCDensitySqDeg,'x','Color',p.Results.cardinalMeridianPlotColors{mm});
    ylim([1,3000]);
    xlim([0.2,70]);
    hold on
    regularSupportPosDegVisual=1e-2:0.01:70;
    loglog(regularSupportPosDegVisual,RGCDensityFit(regularSupportPosDegVisual),'-','Color',p.Results.cardinalMeridianPlotColors{mm});
    xlabel('log10 Eccentricity [deg retina]');
    ylabel('log10 RGC density [counts / deg retina ^2]');
end
legend(num2str(p.Results.cardinalMeridianAngles),'Location','southwest')

subplot(1,2,2)
interpolarMeridianAngles=[45 135 225 315];
for mm=1:4
    [RGCDensityFit] = getSplineFitToRGCDensitySqDegVisual(interpolarMeridianAngles(mm));
    loglog(regularSupportPosDegVisual,RGCDensityFit(regularSupportPosDegVisual),'-','Color',p.Results.cardinalMeridianPlotColors{mm});
    xlim([0.2,70]);
    ylim([1,3000]);
    xlabel('log10 Eccentricity [deg retina]');
    ylabel('log10 RGC density [counts / deg retina ^2]');
    hold on
end
legend(num2str(interpolarMeridianAngles),'Location','southwest')
if p.Results.saveVectorPlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,'rgcDensitySplineModel.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end


%% Plot the spline cone and rgc fits for one meridian in linear space
mm=3;
figHandle = figure();
subplot(2,1,1)
[coneDensitySqDegVisual, coneNativeSupportPosDegVisual] = loadRawConeDensityByEccen(p.Results.cardinalMeridianAngles(mm));
dispConeNativeSupportPosDeg=coneNativeSupportPosDegVisual;
[fitConeDensitySqDegVisual] = getSplineFitToConeDensitySqDegVisual(p.Results.cardinalMeridianAngles(mm));
plot(dispConeNativeSupportPosDeg,coneDensitySqDegVisual,'x','Color',p.Results.cardinalMeridianPlotColors{mm});
xlim([0,30]);
ylim([0,3e4]);
hold on
regularSupportPosDegVisual=0:0.01:70;
plot(regularSupportPosDegVisual,fitConeDensitySqDegVisual(regularSupportPosDegVisual),'-','Color',p.Results.cardinalMeridianPlotColors{mm});
xlabel('Eccentricity [deg retina]');
ylabel('Cone density [counts / deg retina ^2]');
subplot(2,1,2)
[RGCDensitySqDeg, RGCNativeSupportPosDeg] = loadRawRGCDensityByEccen(p.Results.cardinalMeridianAngles(mm));
dispRGCNativeSupportPosDeg=RGCNativeSupportPosDeg;
[RGCDensityFit] = getSplineFitToRGCDensitySqDegVisual(p.Results.cardinalMeridianAngles(mm));
plot(dispRGCNativeSupportPosDeg,RGCDensitySqDeg,'x','Color',p.Results.cardinalMeridianPlotColors{mm});
xlim([0,30]);
ylim([0,2500]);
hold on
regularSupportPosDegVisual=1e-2:0.01:70;
loglog(regularSupportPosDegVisual,RGCDensityFit(regularSupportPosDegVisual),'-','Color',p.Results.cardinalMeridianPlotColors{mm});
xlabel('Eccentricity [deg]');
ylabel('RGC density [counts / deg retina ^2]');
if p.Results.saveVectorPlots
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
[fitConeDensitySqDegVisual] = getSplineFitToConeDensitySqDegVisual(p.Results.cardinalMeridianAngles(mm));
regularSupportPosDegVisual=0:0.01:70;
[ mRFDensitySqDeg ] = transformConeToMidgetRFDensity( fitConeDensitySqDegVisual(regularSupportPosDegVisual), 'linkingFuncParams', fitParamsByMeridian(mm,1:2) );
plot(regularSupportPosDegVisual,mRFDensitySqDeg,'-','Color',p.Results.cardinalMeridianPlotColors{mm});
xlim([0,30]);
ylim([0,3e4]);
xlabel('Eccentricity [deg retina]');
ylabel('mRF density [counts / deg retina ^2]');

subplot(2,1,2)
[RGCDensityFit] = getSplineFitToRGCDensitySqDegVisual(p.Results.cardinalMeridianAngles(mm));
[ mRGCDensitySqDeg ] = transformRGCToMidgetRGCDensityDacey( regularSupportPosDegVisual, RGCDensityFit(regularSupportPosDegVisual), 'linkingFuncParams', fitParamsByMeridian(mm,3:end) );
plot(regularSupportPosDegVisual,mRGCDensitySqDeg,'-','Color',p.Results.cardinalMeridianPlotColors{mm});
xlim([0,30]);
ylim([0,2500]);
xlabel('Eccentricity [deg retina]');
ylabel('mRGC density [counts / deg retina^2]');
if p.Results.saveVectorPlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,'mRGCmRF_linearOneMeridian.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end


%% Plot the mRGC fraction for the cardinal meridians
figHandle = figure();
for mm = 1:4
    
    % Get the spline fit function to the Curcio & Allen RGC density
    fitRGCDensitySqDegVisual = getSplineFitToRGCDensitySqDegVisual( p.Results.cardinalMeridianAngles(mm) );

    % Define a high-resolution eccentricity support
    fitSupportPosDeg = 0:0.01:40;
    
    % Plot the Drasdo and Dacey midget fraction
    subplot(1,2,1);
    f0 = 0.8928; rm = 41.03; % Drasdo's values
    midgetFraction_Drasdo = calcDrasdoMidgetFractionByVisualEccen(fitSupportPosDeg,f0,rm);
    plot(fitSupportPosDeg,midgetFraction_Drasdo,'-k');
    hold on
    midgetFraction_Dacey = calcDaceyMidgetFractionByEccenDegVisual(fitSupportPosDeg);
    plot(fitSupportPosDeg,midgetFraction_Dacey,'-k');
    xlabel('eccentricity visual deg (at the photoreceptors - Drasdo, at the retina - dacey)');
    ylabel('midget fraction');
    ylim([0 1]);
    xlim([0 40]);
    title('Dacey and Drasdo midget fraction');
    pbaspect([2 1 1]);
    
    % Get our midget fraction
    [ ~, midgetFraction_ours ] = transformRGCToMidgetRGCDensityDacey( fitSupportPosDeg, fitRGCDensitySqDegVisual(fitSupportPosDeg), 'linkingFuncParams', fitParamsByMeridian(mm,3:end) );

    % Nan the optic disc
    midgetFraction_ours = nanOpticDiscPoints(midgetFraction_ours, fitSupportPosDeg, p.Results.cardinalMeridianAngles(mm));

    % Plot it
    subplot(1,2,2);
    plot(fitSupportPosDeg,midgetFraction_ours','-','Color',p.Results.cardinalMeridianPlotColors{mm});
    hold on
    ylim([0 1]);
    xlim([0 40]);
    xlabel('eccentricity deg retina');
    ylabel('midget fraction');
    title('Our midget fraction');
    pbaspect([2 1 1]);
    drawnow
end
if p.Results.saveVectorPlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,'mRGCFractionByMeridian.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end

%% Plots related to midget RF density

figHandle = figure();
for mm = 1:4
    
    % Get the spline fit function to the Curcio & Allen RGC density
    fitConeDensitySqDegVisual = getSplineFitToConeDensitySqDegVisual( p.Results.cardinalMeridianAngles(mm) );

    % Define a high-resolution eccentricity support
    fitSupportPosDeg = 0.01:0.01:60;
    
    % Get the vector of cone density
    coneDensitySqDegVisual = fitConeDensitySqDegVisual(fitSupportPosDeg);
    
    % Nan the optic disc
    coneDensitySqDegVisual = nanOpticDiscPoints(coneDensitySqDegVisual, fitSupportPosDeg, p.Results.cardinalMeridianAngles(mm));
    
    % calculate the mRF density using Watson equation 8 
    [ mRFDensitySqDegVisual_watson ] = calcWatsonMidgetRFDensityByEccenDegVisual(fitSupportPosDeg, p.Results.cardinalMeridianAngles(mm));
    
    % Plot the mRF density by eccentricity for Watson
    subplot(2,2,1);
    
    loglog(fitSupportPosDeg,mRFDensitySqDegVisual_watson,'-','Color',p.Results.cardinalMeridianPlotColors{mm});
    ylim([1e0 1e5]);
    xlim([1e-1 1e2]);
    xlabel('log10 eccentricity deg visual');
    ylabel('log10 mRF density / deg visual ^2');
    title('Watson''s mRF density by eccentricity');
    pbaspect([2 1 1]);
    hold on
    
    % Plot the mRF density by eccentricity from our functions
    subplot(2,2,2);
    mRFDensitySqDegVisual_ours = transformConeToMidgetRFDensity(coneDensitySqDegVisual, ...
        'linkingFuncParams',fitParamsByMeridian(mm,1:2));
    
    loglog(fitSupportPosDeg,mRFDensitySqDegVisual_ours,'-','Color',p.Results.cardinalMeridianPlotColors{mm});
    ylim([1e0 1e5]);
    xlim([1e-1 1e2]);
    xlabel('log10 eccentricity visual deg');
    ylabel('log10 mRF density / visual deg2');
    title('Our mRF density by visual eccentricity');
    pbaspect([2 1 1]);
    hold on
    
    % Plot the mRF : cone ratio for Watson
    subplot(2,2,3);
    loglog(fitSupportPosDeg,mRFDensitySqDegVisual_watson./coneDensitySqDegVisual,'-','Color',p.Results.cardinalMeridianPlotColors{mm});
    ylim([1e-2 10]);
    xlim([1e-1 1e2]);
    xlabel('log10 eccentricity visual degrees');
    ylabel('log10 mRF:cone');
    title('Watson''s mRF:cone ratio by visual eccentricity');
    pbaspect([2 1 1]);
    hold on
    
    
    % Plot the mRF : cone ratio for us
    subplot(2,2,4);
    loglog(fitSupportPosDeg,mRFDensitySqDegVisual_ours./coneDensitySqDegVisual,'-','Color',p.Results.cardinalMeridianPlotColors{mm});
    ylim([1e-2 10]);
    xlim([1e-1 1e2]);
    xlabel('log10 eccentricity visual degrees');
    ylabel('log10 mRF:cone');
    title('Our mRF:cone ratio by visual eccentricity');
    hold on
    pbaspect([2 1 1]);
    drawnow
    
end
if p.Results.saveVectorPlots
    fileOutPath = fullfile(p.Results.pathToPlotOutputDirRoot,p.Results.subjectName,'mRFDensityAndRatioByMeridian.pdf');
    saveas(figHandle,fileOutPath)
    close(figHandle);
end

end % function


%% LOCAL FUNCTIONS

function vectorOut = zeroOpticDiscPoints(vectorIn, regularSupportPosDegVisual, polarAngle)
	opticDiscIndices = findOpticDiscPositions(regularSupportPosDegVisual, polarAngle);
    vectorOut = vectorIn;
    vectorOut(opticDiscIndices) = 0;
end