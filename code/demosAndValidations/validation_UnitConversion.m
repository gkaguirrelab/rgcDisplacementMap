function validation_UnitConversion(varargin)
% validation_UnitConversion
%
% Tests if the unit conversion routines are close to invertible
%

%% Parse input and define variables
p = inputParser;

% Optional anaysis params
p.addParameter('cardinalMeridianAngles',[0 90 180 270],@isnumeric);
p.addParameter('errorTolerance',1e-3,@isnumeric);

% parse
p.parse(varargin{:})


%% Clean up
close all

passedAllTestsFlag = true;

%% Test degRetina --> mmRetina --> degVisual and back
for mm = 1:length(p.Results.cardinalMeridianAngles)
    supportPosDegRetinaIn = 0:1:100;
    vectorInTranslation = convert_degRetina_to_mmRetina(supportPosDegRetinaIn);
    vectorInTranslation = convert_mmRetina_to_degVisual(vectorInTranslation, p.Results.cardinalMeridianAngles(mm));
    vectorInTranslation = convert_degVisual_to_mmRetina(vectorInTranslation, p.Results.cardinalMeridianAngles(mm));
    supportPosDegRetinaOut = convert_mmRetina_to_degRetina(vectorInTranslation);
    if max(abs(supportPosDegRetinaIn-supportPosDegRetinaOut)./supportPosDegRetinaIn) > p.Results.errorTolerance
        fprintf('Inversion error in position units greater than a proportion of %d \n',p.Results.errorTolerance);
        passedAllTestsFlag = false;
    end
end


%% Test mmSqRetina to degSqVisual and bacl
for mm = 1:length(p.Results.cardinalMeridianAngles)
    supportPosMMRetina = 0:1:25;
    valuesPer_mmSqRetinaIn = rand(size(supportPosMMRetina));
    factorA = calc_degSqVisual_per_mmSqRetina(supportPosMMRetina, p.Results.cardinalMeridianAngles(mm));
    vectorInTranslation = valuesPer_mmSqRetinaIn ./ factorA;
    factorB = calc_mmSqRetina_per_degSqVisual(convert_mmRetina_to_degVisual(supportPosMMRetina, p.Results.cardinalMeridianAngles(mm)), p.Results.cardinalMeridianAngles(mm));
    valuesPer_mmSqRetinaOut = vectorInTranslation ./ factorB;
    if max(abs(valuesPer_mmSqRetinaIn-valuesPer_mmSqRetinaOut)./valuesPer_mmSqRetinaIn) > p.Results.errorTolerance
        fprintf('Inversion error in area units greater than a proportion of %d \n',p.Results.errorTolerance);
        passedAllTestsFlag = false;
    end
end

if passedAllTestsFlag
    fprintf('Passed all unit conversion tests\n');
end


end % function