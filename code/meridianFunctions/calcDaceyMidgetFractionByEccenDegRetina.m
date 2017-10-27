function midgetFraction = calcDaceyMidgetFractionByEccenDegRetina(supportPosDegRetina, varargin)
% calcDaceyMidgetFractionByEccen - midget fraction as a function of eccentricity
%
% This routine returns, for each of the locations in supportPosDegRetina,
% the fraction of retinal ganglion cells that are midget RGCs. The
% calculation is based upon the values provided in Dacey (1993) J Neurosci.
% A logisitic function is fit to the data from the Dacey paper which then
% allows for interpolation to intermediate values.
%

%% Parse input and define variables
p = inputParser;

% required input
p.addRequired('supportPosDegRetina',@isnumeric);

% Optional anaysis params
p.addParameter('minRatio',0.45,@isnumeric);
p.addParameter('maxRatio',0.95,@isnumeric);
p.addParameter('logitFitStartPoint',[5, 20],@isnumeric);

% parse
p.parse(supportPosDegRetina, varargin{:})


% Define a four-parameter logistic function that will be used to fit the
% modeled relationship. Two of the parameters (max and min asymptote) are
% locked by the passed parameter
logisticFunc = fittype( @(slope,inflect,minRatio,maxRatio,x) minRatio+(maxRatio-minRatio)./(1+sign(x./inflect).*abs((x./inflect).^slope)), ...
    'independent','x','dependent','y','problem',{'minRatio','maxRatio'});

% Data taken from Figure 19B of:
%	Dacey 1993 J Neurosci,
%   The Mosaic of Midget Ganglion Cells in the Human Retina
% digitzed using the web plot digitizer tool
daceyDataSupportPosMmRetina = [0.98174998	2.02854114	3.015966729	3.978015168	4.956959961	6.011644785	6.96599527	7.975927587	9.039745576	10.0028378	10.99287287	11.98068988	13.01939167	14.00851341	14.99998369];
daceyMidgetFraction = [91.54203702	97.96884938	95.9768409	94.73766615	82.97154041	68.49580038	58.38473457	52.33172959	48.381962	48.34575552	49.36116774	47.82027236	44.92408057	44.88689554	47.55638914];

% Convert the Dacey support vector from mm to degrees
daceyDataSupportPosDegRetina = convert_mmRetina_to_degRetina(daceyDataSupportPosMmRetina);

% Convert the Dacey percent values to proportion
daceyMidgetFraction = daceyMidgetFraction / 100;

% Fit a logistic function
logisticFit = fit(daceyDataSupportPosDegRetina',daceyMidgetFraction',logisticFunc, ...
    'problem',{p.Results.minRatio, p.Results.maxRatio}, ...
    'StartPoint',p.Results.logitFitStartPoint, ...
    'Lower',[0,0],'Upper',[100,100] );

midgetFraction = logisticFit(supportPosDegRetina);

end