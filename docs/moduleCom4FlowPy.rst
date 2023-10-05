com4FlowPy: Flow-Py
==================

.. Note::
  THIS MODULE IS CURRENTLY UNDER HEAVY DEVELOPMENT! 

  The code is not automatically tested in any way and not included in the code coverage!

  It also does not adhere to the AvaFrame coding and naming conventions (yet)...

  Use at your own risk and if you want to contribute to the modules improvement, you are very welcome to do so!


Flow-Py is an open source tool to compute gravitational mass flows (GMF) run out and intensity. 
The main objective is to compute the spatial extent of GMF, which consists of the starting, 
transit and runout zones of GMF on the surface of a three dimensional terrain. The resulting 
run out is mainly dependent on the terrain and the location of the starting/release point.
No time-dependent equations are solved in the model. com4FlowPy uses existing statistical-data-based
approaches for solving the routing and stopping of GMF. 

The tool has been designed to be computationally light, allowing the application on a 
regional scale including a large number of GMF paths. The 
implementation allows users to address specific GMF research
questions by keeping the parameterization flexible and the ability to include
custom model extensions and add-ons.

Running the Code
----------------

**You are required to install rasterio and gdal separately, since they are not included 
in the general AvaFrame requirements.**

If you have trouble installing GDAL or rasterio on Windows use these links to
get the required version directly from their website, first install *GDAL* and then *rasterio*.

*GDAL*: https://www.lfd.uci.edu/~gohlke/pythonlibs/#gdal

*rasterio*: https://www.lfd.uci.edu/~gohlke/pythonlibs/#rasterio

Once the required libraries are installed the model runs via the ``runCom4FlowPy.py`` script. 

Configuration
----------------

The configuration can be found in ``com4FlowPyCfg.ini``

- alpha_angle (controls the run out angle induced stopping and routing)
- exponent (controls concentration of routing flux and therefore the lateral spread)
- working directory path
- path to DEM raster (.tiff or .asc)
- path to release raster (.tiff or .asc)  
- (Optional) flux threshold (positive number) flux_threshold=xx (limits spreading with the exponent)
- (Optional) Max Z<sup>&delta;</sup> (positive number) max_z_delta=xx (max kinetic energy height, turbulent friction)

Input Files
-----------

All raster files (DEM, release, ...) must be in the .asc or .tif format.

All rasters need the same resolution (normal sizes are e.g. 5x5 or 10x10 meters).

All Layers need the same spatial extend, with no data values < 0 (standard no data values = -9999).

The locations identified as release areas need values > 0. (see release.tif in examples)

Output
------

All outputs are in the .tiff raster format in the same resolution and extent as the input raster layers.

- z_delta: the maximum z_delta of all paths for every raster cell (geometric measure of process magnitude, can be associated to kinetic energy/velocity)
- Flux: The maximum routing flux of all paths for every raster cell
- sum_z_delta: z_delta summed up over all paths on every raster cell
- Cell_Counts: number of paths that route flux through a raster cell
- Flow Path Travel Angle, FP_TA: the gamma angle along the flow path
- Straight Line Travel Angle, SL_TA: Saves the gamma angle, while the distances are calculated via a straight line from the release cell to the current cell


