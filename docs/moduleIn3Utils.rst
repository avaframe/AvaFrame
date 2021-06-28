######################################
in3Utils: Various Utilities Modules
######################################


geoTrans
============
The :py:mod:`in3Utils.geoTrans` module provides useful functions to operate transformations,
comparison or interpolation on raster, lines, points...
Further information about the available functions in :py:mod:`in3Utils.geoTrans`



Generate Topography
===================

The :py:mod:`in3Utils.generateTopo` module provides functions to generate DEM files
for idealized/generic topographies that can be used for snow avalanche simulations.

This module can generate the following topographies:

* flat plane (FP)
* inclined plane of constant slope (IP)
* parabola - parabolic slope transitioning into a flat foreland (PF)
* hockey stick - inclined plane of constant slope with a smooth transition to a flat foreland (HS)
* bowl-shaped topography (BL)
* helix-shaped topography (HX)
* pyramid-shaped topography (PY)

Extra features can be added to the above topographies:

	* a channel can be introduced (then set channel=True).
		This channel can also be set to widen towards the top and bottom of the channel (then set narrowing=True).
		The channel can be added by either adding a 'channel layer'
		(max thickness=channel depth) on top of the topography (topoAdd=True)
		or by cutting them into the original topography (topoAdd=False).
	* In case of the parabola topography, a Dam can be added by setting dam=True.

Detailed information on the individual functions used to create the topographies can be found in:

:py:mod:`in3Utils.generateTopo`

Input
------

* in ``generateTopoCfg.ini`` all required input parameters are listed (includes default values for all parameters).

Output
-------

* 3D surface plot of the generated topography as .png file and possibly shown to screen (see flags)
* .asc file of DEM data


To run
------

* copy ``generateTopoCfg.ini`` to ``local_generateTopoCfg.ini`` and set the desired parameter values (if not, the default values are used)
* in ``Avaframe/`` run::

	python3 runGenerateTopo.py


Theory
------

Topographies are generated using inclined and flat planes, parabolas, spheres and circles.
Channels are introduced as half-sphere shaped features with smooth transition
from the no channel area to the channel using cumulative distribution functions.

Configuration parameters
------------------------

In the case of the pyramid-shaped topography, the domain extent is defined by the max elevation (z0 - elevation of the appex point) and the
slope of the pyramid facets (meanAlpha)

**Domain parameters:**

:dx: DEM spatial resolution [m]
:xEnd: total horizontal extent of the domain [m]
:yEnd: total vertical extent of the domain [m]


**Topography parameters:**

:flens: 	distance to the point where the slope transitions into a flat plane [m]
:meanAlpha: 	slope angle from the max. elevation to the start of the flat plane [°] - or slope of the inclined plane [°]
:C: 		total fall height [m]
:rBowl:	  bowl radius [m]
:rHelix:	  radius for helix [m]
:z0: 		max elevation [m]
:zElev:	  elevation of the flat plane [m]
:rCirc: 	radius of the smoothing circle [m]
:demType: topography types (FP, IP, PF, HS, BL, HX, PY - explanation given in the introductory description)
:flatx: extent of the flat foreland for the pyramid in x
:flaty: extent of the flat foreland for the pyramid in y
:phi: rotation angle for the pyramid


**Flags for channels and plotting:**

:channel: True - introduce channel; False - no channel
:narrowing: True - channel is wide at start and end and narrow in the middle part; False - channel is uniform
:topoAdd: True - add channel layer; False: cut channel into original topography;
:flagRot: True - rotate pyramid along z-axis


**Channel parameters:**

:cRadius:  standard channel radius
:cInit: start and end half width of channel that is narrower in the middle part
:cff: standard deviation sigma
:cMustart: mean mu - represents upper part of the channel
:cMuend: mean mu - represents lower part of the channel


Get Release Area
===================

Generate a release area for a topography created with ``generateTopo`` , this function is available for the following topographies:

* flat plane (FP)
* inclined plane (IP)
* parabola (PF)
* hockey stick (HS)

The release areas are defined as rectangular features build by four corner points, which are based on the following conditions:

* prescribed vertical stretch in meters (difference in altitude)
* prescribed volume of the release area
* lower margin is located where the slope angle falls below 30°
* if slope does not fall below 30 °, the upper margin is located xStart away from the upper margin of the DEM

Detailed information on the individual functions used to create the topographies can be found in :py:mod:`in3Utils.getReleaseArea`


Input
------

* in ``getReleaseAreaCfg.ini`` and ``generateTopoCfg.ini`` all the required input parameters are listed (does include default values for all parameters)

Output
------

* release area as shapefile, .nxyz and .txt file
* if showplot flag is True, plot of release area on domain extent

To run
------

Following these steps, you can generate an avalanche test case including a DEM and a simple realease area.

* copy ``generateTopoCfg`` and  ``getReleaseAreaCfg`` to ``local_generateTopoCfg.ini``
  and ``local_getReleaseAreaCfg.ini`` and set desired parameter values (if not, the default values are used)
* in ``avaframe`` run::

	python3 runGenProjTopoRelease.py

**Parameters:**

:hr: release area vertical stretch [m]
:vol: volume of snow in release area [m3]
:dh: release snow thickness [m]
:xStart: upper margin of the release area distance in x from origin [m]
:lenP: number of release area polygon points
:outputtxt: True - copy the output to txt file
:xExtent: horizontal extent of release area for flat plane
:alphaStop: slope angle that defines lower margin of release area
:relNo: number of release area for name
:relName: name of release area feature in shapefile


Initialize Project
=====================

:py:mod:`in3Utils.InitializeProject` provides functions to initialize a project, create the required directory structure and delete specified
files or directories.

The main function *initializeFolderStruct*, creates the folder structure required to perform avalanche simulations:
::

		NameOfAvalanche/
			Inputs/
				ENT/		- entrainment areas
				LINES/		- avalanche paths
				POINTS/		- split points
				REL/		- release area scenario
				RES/		- resistance areas
				SECREL/ - secondary release areas
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

* copy ``avaframeCfg.ini`` to ``local_avaframeCfg.ini`` and set your desired avalanche directory name
* in ``avaframe`` run:

			python3 runInitializeProject.py


..
   Logging
   ========

   .. include:: moduleLogging.rst

   Configuration
   ==============

   .. include:: moduleConfiguration.rst



fileHandlerUtils
=================

:py:mod:`in3Utils.fileHandlerUtils` gathers useful functions to create directories, read log files,
extract information from logs, fetch and export data and fetch simulation info into a dataFrame
that can be used within other functions. Details on these functions can be found in:

:py:mod:`in3Utils.fileHandlerUtils`
