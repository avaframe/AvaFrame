##################################
in3Utils: Module Project Utilities
##################################


geoTrans
===================
The ``geoTrans.py`` module gathers useful functions to operate transformations on raster, lines, points...

Functions
-----------

**Projection on Raster:**

``Points = projectOnRaster(dem, Points)`` takes a "dem" dictionary and "Points" dictionary
(wich can be a single point or a line...) in input and returns the "Points" dictionary with
an extra "z" argument representing the "z" coordinate of the (x,y) point on the dem.


``Points = projectOnRasterVect(dem, Points, interp = 'bilinear')`` Does the same as the previous
function but also operates on 2D arrays. All the calculation are vectorized to avoid loops.
Two interpolation methods are available, 'nearest' or 'bilinear'


**Prepare Line:**

``AvaProfile, projSplitPoint = prepareLine(dem, AvaPath, distance=10, Point=None)`` takes a "dem" dictionary,
a "AvaPath" dictionary (x, y coordinates), a re-sampling distance and a "Point" dictionary in input and returns
the "AvaProfile" dictionary corresponding to the "AvaPath" dictionary. That is to say the "line" dictionary re-sampled
according to distance with the corresponding "z" argument representing the "z" coordinate of the re-sampled (x,y)
point on the dem and the curvilinear coordinate "s" along the line (the first point of the line has a s=0).
It also returns the projection of the Point on the AvaProfile if this one was supplied in input.

**Project on Profile:**

``projSplitPoint = findSplitPoint(AvaProfile, splitPoint)`` takes a "AvaProfile" dictionary
and a "splitPoint" dictionary in input and returns the "projSplitPoint" dictionary which is the projection of
"splitPoint" on the "AvaProfile".


**Check Profile:**

``projSplitPoint, AvaProfile = checkProfile(AvaProfile, projSplitPoint=None)`` takes a "AvaProfile" dictionary
and a "projSplitPoint" dictionary in input and check if the Profile goes from top to bottom,
reverts it if necessary and returns the correct "projSplitPoint" and "AvaProfile" dictionaries.

**Prepare inputs for find angle in profile:**

``angle, tmp, deltaInd =prepareAngleProfile(beta, AvaProfile)`` takes a angle value in degres and
an Avalanche profile in input, computes the angle of the profile and returns this ``angle``, the list
of indexes ``tmp`` where the angle is under the input angle value and ``deltaInd`` the number of consecutive
indexes required.

**Find angle in profile:**

``idsAnglePoint =findAngleProfile(tmp, deltaInd)`` takes the outputs of ``prepareAngleProfile`` as inputs
and returns the index of the desired angle as output.

**Bresenham Algorithm:**

``z = bresenfindCellsCrossedByLineBresenhamham(x0, y0, x1, y1, cs)`` takes a two (x,y) points and a cell size in input and returns
the z = (x,y) list of the cells hit by the line between the two input points.


**Path to domain:**

``rasterTransfo = path2domain(xyPath, rasterTransfo)`` takes the (x,y) coordinates of a polyline,
a domain width and a cell size (in rasterTransfo) in input and returns the domain of width w along the polyline.

**Polygon to mask:**

``mask = poly2mask_simple(ydep, xdep, ncols, nrows)`` takes the (x,y) coordinates
of a polygon and a rater size in input and returns the raster mask corresponding to the polygon.

**In polygon:**

``IN = inpolygon(X, Y, xv, yv)`` takes the (X, Y) coordinates of points and xv, yv foot print of a
polygon on a raster in input and returns the raster mask corresponding to the polygon.



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
------

* in ``generateTopoCfg.ini`` all required input parameters are listed (does include default values for all parameters)

Output
-------

* 3D surface plot of generated topography as .png file and possibly shown to screen (see flags)
* .asc file of DEM data


To run
------

* copy ``generateTopoCfg.ini`` to ``local_generateTopoCfg.ini`` and set desired parameter values (if not, the default values are used)
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
:xEnd: total horizontal extent of the domain [m]
:yEnd:	 total horizontal extent of the domain [m]


**Topography parameters:**

:f_lens: 	distance to point where slope transitions into flat plane [m]
:mean_alpha: 	slope angle from max. elevation to start flat plane [°] - or slope of inclined plane [°]
:C: 		total fall height [m]
:rBowl:	  bowl radius [m]
:rHelix:	  radius for helix [m]
:z0: 		max elevation [m]
:zElev:	  elevation of flat plane [m]
:rCirc: 	radius of smoothing circle [m]
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


Get Release Area
===================

Generate a release area for a topography created with ``generateTopo`` , this function is available for the following topographies:

* flat plane (FP)
* inclined plane (IP)
* hockey stick (HS)
* hockey stick smooth (HS2)

The release areas are defined as rectangular features build by four corner points, which are based on the following conditions:

* prescribed vertical stretch of 200 m (difference in altitude)
* prescribed volume of the release area
* lower margin is located where the slope angle falls below 30°
* if slope does not fall below 30 °, upper margin is located xStart away from upper margin of the DEM

Input
------

* in ``getReleaseAreaCfg.ini`` and ``generateTopoCfg.ini`` all required input parameters are listed (does include default values for all parameters)

Output
------

* release area as shapefile, .nxyz and .txt file
* if showplot flag is True, plot of release area on domain extent

To run
------

Following this steps, you can generate an avalanche test case including a DEM and a simple realease area.

* copy ``generateTopoCfg`` and  ``getReleaseAreaCfg`` to ``local_generateTopoCfg.ini``
  and ``local_getReleaseAreaCfg.ini`` and set desired parameter values (if not, the default values are used)
* in ``avaframe`` run::

	python3 runGenProjTopoRelease.py

**Parameters:**

:hr: release area vertical stretch [m]
:vol: volume of snow in release area [m3]
:dh: release snow thickness [m]
:xStart: upper margin of release area distance in x from origin [m]
:lenP: number of release area polygon points
:showPlot: True - show plot of release area
:outputtxt: True - copy the output to txt file
:xExtent: horizontal extent of release area for flat plane
:alphaStop: slope angle that defines lower margin of release area
:relNo: number of release area for name
:relName: name of release area feature in shapefile


Initialize Project
=====================

This function creates the folder structure required to perform avalanche simulations:
::

		NameOfAvalanche/
			Inputs/
				REL/		- release area scenario
				RES/		- resistance areas
				ENT/		- entrainment areas
				POINTS/		- split points
				LINES/		- avalanche paths
				.asc		- DEM
			Outputs/
			Work/


Input
-------

* path to NameOfAvalanche

This path is specified in the configuration file ``avaframeCfg.ini`` with the parameter *avalancheDir*.

Output
-------

* NameOfAvalanche directory


To run
-------

*  copy ``avaframeCfg.ini`` to ``local_avaframeCfg.ini`` and set your desired avalanche directory name
* in ``avaframe`` run:

			python3 runInitializeProject.py


Logging
========

.. include:: moduleLogging.rst

Configuration
==============

.. include:: moduleConfiguration.rst
