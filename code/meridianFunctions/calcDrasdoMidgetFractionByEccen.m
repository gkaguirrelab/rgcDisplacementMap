function midgetFraction = calcDrasdoMidgetFractionByEccen(supportPosDeg,f0,rm)

% The equation is taken from Drasdo et al 2007.
midgetFraction = f0.*(1+(supportPosDeg./rm)).^-1;

end