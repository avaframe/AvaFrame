##################################
in3Utils: Module Project Utilities
##################################



Generate Topography
===================

Generate DEM files for idealised/generic topographies that can be used as required DEM input data for snow avalanche simulations.

The following topography types can be chosen:

* flat plane (FP)
* inclined plane of constant slope (IP)
* hockey stick - using a parabola to describe the sloping plane that then transitions into a flat foreland (HS)
* hockey stick smooth - using an linearly sloping plane with smoothed transition to the flat foreland (HS2)
* bowl-shaped topography (BL)
* helix-shaped topography (HX)

On top of these topographies, channels can be introduced (then set flag_channel=True), and these channels can also be set to become narrower along the channel and wider at the end of the channel (then set flag_narrowing=True).
There is the option to introduce these channels by either cutting them into the original topography (flag_topoconst=True) or by adding a 'channel layer' (max thickness=channel depth) on top of the topography.

Input
-----

* in ``generateTopoCfg.ini`` all required input parameters are listed (does include default values for all parameters)

Outputs
-------

* 3D surface plot of generated topography as .png file and possibly shown to screen (see flags)
* .asc file of DEM data


To run
------

* copy ``generateTopoCfg.ini`` to ``local_generateTopoCfg.ini`` (and set desired parameter values, if not, the default values are used)
* in ``Avaframe/`` run::

	python3 runGenerateTopo.py


Theory
------

Topographies are generated using inclined and flat planes, parabolas, spheres and circles.
Channels are introduced as half-sphere shaped features with smooth transition from no channel to channel using cumulative distribution functions.

Configuration parameters
------------------------

**Domain parameters:**

:dx: DEM spatial resolution [m]
:x_end: total horizontal extent of the domain [m]
:y_end:	 total horizontal extent of the domain [m]


**Topography parameters:**

:f_lens: 	distance to point where slope transitions into flat plane [m]
:mean_alpha: 	slope angle from max. elevation to start flat plane [°] - or slope of inclined plane [°]
:C: 		total fall height [m]
:r_bowl:	  bowl radius [m]
:r_helix:	  radius for helix [m]
:z0: 		max elevation [m]
:z_elev:	  elevation of flat plane [m]
:r_circ: 	radius of smoothing circle [m]
:DEM_type: topography types (FP, IP, HS, HS2, BL, HX - explanation given in the introductory description**


**Flags for channels and plotting:**

:channel: True - introduce channel; False - no channel
:narrowing: True - channel is wide at start and end and narrow in the middle part; False - channel is uniform
:topoconst: True - cut channel into original topography; False - add channel layer
:showplot: True - show plot of generated DEM

**Channel parameters:**

:c_radius:  standard channel radius
:c_init: start and end half width of channel that is narrower in the middle part
:c_ff: standard deviation sigma
:c_mustart: mean mu - represents upper part of the channel
:c_muend: mean mu - represents lower part of the channel


Read / Write / Compare raster from or to ASCII files
=========================================================

``ascUtilis.py`` is a module created to handle raster ASCII files. It contains different functions
to read ASCII files, either just the header or also the raster matrix and write the data to a numpy array or to
compare raster file headers as well as to write a raster to an ASCII file given a header and data.

Functions
------------------------

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
header information as well as the raster data in a numpy array.

**Write ASCII file:**

``writeResultToAsc(header, resultArray, outType=None):`` takes a header and numpy array as inputs and writes the
corresponding raster ASCII file.
