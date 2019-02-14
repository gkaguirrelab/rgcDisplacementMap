function vectorOut = zeroOpticDiscPoints(vectorIn, regularSupportPosDegVisual, polarAngle)
	opticDiscIndices = findOpticDiscPositions(regularSupportPosDegVisual, polarAngle);
    vectorOut = vectorIn;
    vectorOut(opticDiscIndices) = 0;
end
