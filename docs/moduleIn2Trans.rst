##################################
in2Trans: Transformation Utilities
##################################

Working with ASCII files
==========================

:py:mod:`in2Trans.ascUtils` is a module created to handle raster ASCII files. It
contains different functions to read ASCII files and write the data to a numpy
array, to compare raster file headers or to write a raster to an ASCII file.
A description of the functions is available in
:py:mod:`in2Trans.ascUtils`.


Working with shapefiles
=============================

:py:mod:`in2Trans.shpConversion` is a module for handling shapefiles.
It contains different functions to read shapefiles and convert them to a python dictionary.
It also provides functions to extract or remove a feature from the shapefile dictionary.


Reading shapefiles
--------------------

Shapefiles are converted to a python dictionary. The dictionary has information
about the number of features as well as the coordinates of the points.
The output dictionary SHPdata looks like this::

		SHPdata['Name'] = ['nameFeature1', 'nameFeature2', 'nameFeature3']
		SHPdata['x'] = [xCoordsssFeature1, xCoordsssFeature2, xCoordsssFeature3]
		SHPdata['y'] = [yCoordsssFeature1, yCoordsssFeature2, yCoordsssFeature3]
		SHPdata['z'] = [zerosss, zerosss, zerosss]
		SHPdata['Start'] = [indexStartFeature1, indexStartFeature2, indexStartFeature3]
		SHPdata['Length'] = [lenghtFeature1, lenghtFeature2, lenghtFeature3]
		SHPdata['sks'] = 'ProjectionInformation'

A description of the functions is available in :py:mod:`in2Trans.shpConversion`
