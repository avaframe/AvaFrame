##################################
in2Trans: Transformation Utilities
##################################

Working with ASCII files
=========================================================

``ascUtilis.py`` is a module created to handle raster ASCII files. It contains different functions
to read ASCII files, either just the header or also the raster matrix and write the data to a numpy array or to
compare raster file headers as well as to write a raster to an ASCII file given a header and data.

.. _ascii:

Functions
------------

**Header Class:**

Header information are read and stored as object from ``class cASCheader``:
::

		header = cASCheader()
		header.nrows = None
		header.ncols = None
		header.cellsize = None
		header.xllcorner = None
		header.xllcenter = None
		header.yllcorner = None
		header.yllcenter = None
		header.noDataValue = None

**Read ASCII header:**

``header = readASCheader(fname)`` takes a .asc file name as input and returns the header information.

**Compare headers:**

``isEqualASCheader(headerA, headerB)`` takes two headers as input and checks if they are equal.

**Read ASCII file to numpy array:**

``header, data = readASCdata2numpyArray(fName, headerFile=None):`` takes a .asc file name as input and returns the
raster data in a numpy array.


**Read ASCII file to dictionary:**

``rasterHeaderData = readRaster(fname):`` takes a .asc file name as input uses readASCdata2numpyArray and returns the
header information as well as the raster data in a numpy array in a dictionary.
::

		rasterHeaderData['header'] = header
		rasterHeaderData['rasterData'] = rasterData


**Write ASCII file:**

``writeResultToAsc(header, resultArray, outType=None):`` takes a header and numpy array as inputs and writes the
corresponding raster ASCII file.


Working with shape files
=============================

``shpConversion.py`` is a module created to handle shape files. It contains different functions
to read shape files to numpy arrays, either lines or points


.. _shape:

Functions
-------------

**Read shape file:**

``SHPdata = SHP2Array(fname, defname=None)`` takes a .shp file name as input (and eventualy a default name for the layer)
and returns a SHPdata dictionnary containing the layer information (can be multiple points or lines):
::

		SHPdata['Name'] = 'list of paths names'
		SHPdata['x'] = 'np array of the x coords of points in paths'
		SHPdata['y'] = 'np array of the y coords of points in paths'
		SHPdata['z'] = 'np array of the z coords of points in paths'
		SHPdata['Start'] = 'list of starting index of each Line in 'x''
		SHPdata['Length'] = 'list of length of each Line in 'x''

**Read shape file as Lines:**

``Line = readLine(fname, defname, header)`` takes a .shp file name as input,  a default name for the layer and a DEM header
reads the shape file, checks that the Lines lay on the DEM and returns the SHPdata dictionnary containing the Lines information.


**Read shape file as Points:**

``Points = readPoints(fname, header)`` takes a .shp file name as input,  a default name for the layer and a DEM header
reads the shape file, checks that the Lines lay on the DEM and returns the SHPdata dictionnary containing the Points information.
