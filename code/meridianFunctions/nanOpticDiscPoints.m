function vectorOut = nanOpticDiscPoints(vectorIn, regularSupportPosDegVisual, polarAngle)
	opticDiscIndices = findOpticDiscPositions(regularSupportPosDegVisual, polarAngle);
    vectorOut = vectorIn;
    vectorOut(opticDiscIndices) = nan;
end
