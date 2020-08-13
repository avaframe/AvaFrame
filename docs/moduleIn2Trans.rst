##################################
in2Trans: Module Project Utilities
##################################



geoTrans
===================



Reading shape files
=============================

``shpConversion.py`` is a module created to handle shape files. It contains different functions
to read shape files to numpy arrays, either lines or points

Functions
------------------------

**Read shape file:**

``SHPdata = SHP2Array(fname, defname=None)`` takes a .shp file name as input (and eventualy a default name for the layer)
and returns a SHPdata dictionnary containing the layer information (can be multiple points or lines):
::

		SHPdata['Name'] = list of paths names
		SHPdata['x'] = np array of the x coords of points in paths
		SHPdata['y'] = np array of the y coords of points in paths
		SHPdata['z'] = np array of the z coords of points in paths
		SHPdata['Start'] = list of starting index of each Line in 'x'
		SHPdata['Length'] = list of length of each Line in 'x'

**Read shape file as Lines:**

``Line = readLine(fname, defname, header)`` takes a .shp file name as input,  a default name for the layer and a DEM header
reads the shape file, checks that the Lines lay on the DEM and returns the SHPdata dictionnary containing the Lines information.


**Read shape file as Points:**

``Points = readPoints(fname, header)`` takes a .shp file name as input,  a default name for the layer and a DEM header
reads the shape file, checks that the Lines lay on the DEM and returns the SHPdata dictionnary containing the Points information.
