function midgetFraction = calcDaceyMidgetFractionByEccen(supportPosDeg)

% Define a four-parameter logistic function that will be used to fit the
% modeled relationship. Two of the parameters (max and min asymptote) are
% locked by the passed parameter
logisticFunc = fittype( @(slope,inflect,minRatio,maxRatio,x) minRatio+(maxRatio-minRatio)./(1+sign(x./inflect).*abs((x./inflect).^slope)), ...
    'independent','x','dependent','y','problem',{'minRatio','maxRatio'});

% Data taken from Figure 19B of:
%	Dacey 1993 J Neurosci,
%   The Mosaic of Midget Ganglion Cells in the Human Retina
% digitzed using the web plot digitizer tool
daceyDataSupportPosMM = [0.98174998	2.02854114	3.015966729	3.978015168	4.956959961	6.011644785	6.96599527	7.975927587	9.039745576	10.0028378	10.99287287	11.98068988	13.01939167	14.00851341	14.99998369];
daceyMidgetFraction = [91.54203702	97.96884938	95.9768409	94.73766615	82.97154041	68.49580038	58.38473457	52.33172959	48.381962	48.34575552	49.36116774	47.82027236	44.92408057	44.88689554	47.55638914];

% Convert the Dacey support vector to degrees
daceyDataSupportPosDeg = convert_mm_to_deg(daceyDataSupportPosMM);

% Convert the Dacey % values to proportion
daceyMidgetFraction = daceyMidgetFraction / 100;

% Fit a logistic function
    blah = fit(daceyDataSupportPosDeg',daceyMidgetFraction',logisticFunc, ...
        'problem',{0.45, 0.95}, ...
        'StartPoint',[5, 20], ...
        'Lower',[-100,0],'Upper',[100,100] );

midgetFraction = blah(supportPosDeg);
    
end