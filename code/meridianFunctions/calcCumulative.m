
function countsPerRing = calcCumulative(regularSupportPosDeg, densityFunction)

extendedRegularSupportPosDeg = [regularSupportPosDeg (regularSupportPosDeg(end)+diff(regularSupportPosDeg(1:2)))];
ringArea = diff(extendedRegularSupportPosDeg.^2 * pi);
countsPerRing = cumsum(densityFunction.*ringArea);

end