in3Utils: Module Project Utilities - generateTopo.py
====================================================

Generate DEM files for idealised/generic topographies that can be used as required DEM input data for snow avalanche simulations.

The following topography types can be chosen:

* flat plane
* inclined plane of constant slope
* hockey stick - using a parabola to describe the sloping plane that then transitions into a flat foreland
* hockey stick smooth - using an linearly sloping plane with smoothed transition to the flat foreland
* bowl-shaped topography
* helix-shaped topography

On top of these topographies, channels can be introduced (then set flag_channel=1), and these channels can also be set to become narrower along the channel and wider at the end of the channel (then set flag_narrowing=1).
There is the option to introduce these channels by either cutting them into the original topography or by adding a 'channel layer' (max thickness=channel depth) on top of the topography.

Input
-----

* in generateTopoCfg.ini all required input parameters are listed (does include default values for all parameters)

Outputs:
--------

* 3D surface plot of generated topography
* .asc file of DEM data

To run:
-------

* copy generateTopoCfg.ini to local_generateTopoCfg.ini (and set desired parameter values, if not, the default values are used)
* in Avaframe/ run:
	python3 runGenerateTopo.py


Theory:
-------
Topographies are generated using inclined and flat planes, parabolas, spheres and circles.
Channels are introduced as half-sphere shaped features with smooth transition from no channel to channel using cumulative distribution functions.

Configuration parameters:
-------------------------

Domain parameters:

* dx:  		DEM spatial resolution [m]
* x_end:	 total horizontal extent of the domain [m]
* y_end:	 total horizontal extent of the domain [m]


Topography parameters:

* f_lens: 	distance to point where slope transitions into flat plane [m]
* mean_alpha: 	slope angle from max. elevation to start flat plane [°] - or slope of inclined plane [°]
* C: 		total fall height [m]
* r_bowl:	  bowl radius [m]
* r_helix:	  radius for helix [m]
* z0: 		max elevation [m]
* z_elev:	  elevation of flat plane [m]
* r_circ: 	radius of smoothing circle [m]
* DEM_type


Flags for channels:

* channel: 1 - introduce channel; 0 - no channel
* narrowing: 1 - channel is wide at start and end and narrow in the middle part; 0 - channel is uniform
* topoconst: 1 - cut channel into original topography; 0 - add channel layer


Channel parameters:

* c_radius:  standard channel radius
* c_init: start and end radius of channel that is narrower in the middle and part
* c_ff: standard deviation sigma
* c_mustart: mean mu - represents upper part of the channel
* c_mustart: mean mu - represents lower part of the channel
