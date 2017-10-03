The Curcio and Allen data contained in this data were retrieved from:

	https://info.cis.uab.edu/curcio/GanglionCellTopography/

Last checked June 28, 2017

Data from this paper:

Curcio, Christine A., and Kimberly A. Allen. "Topography of ganglion cells in human retina." Journal of comparative Neurology 300.1 (1990): 5-25.

The file 'curcioRGCDensityPerSqMm.mat' contains the mean RGC density values (and associated support in mm eccentricity) from these data.

Curcio's measurements were made on the retinal surface in units of mm and counts per mm^2. While the source data provides an eccentricity support vector in units of degrees, this differs a bit from the values produced by our mm --> deg conversion function (which in turn is derived from Watson's paper). So, we load the data in mm / mm^2 and then convert within the routine.