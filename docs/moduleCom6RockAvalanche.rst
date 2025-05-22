com6RockAvalanche: Rock Avalanche
=================================

.. Warning:: This is highly experimental and not tested!

The com6RockAvalanche computational module provides an override setting for com1DFA targeting the simulation of rock
avalanches.

Input
-------

The standard inputs required to perform a simulation run using :py:mod:`com1DFA` 
can be found here: :ref:`moduleCom1DFA:Input`.
However there is one main difference: com6RockAvalanche NEEDS a release thickness raster file. This file has to have
the exact same dimensions as the topography file.
There is a run script to perform a rock avalanche com1DFA run: :py:mod:`runCom6RockAvalanche.py`,
and the configuration settings can be found in ``com6RockAvalanche/com6RockAvalancheCfg.ini``.

To run
------

* first go to ``AvaFrame/avaframe``
* copy ``avaframeCfg.ini`` to ``local_avaframeCfg.ini`` and set your desired avalanche directory name
* create an avalanche directory with required input files - for this task you can use :ref:`moduleIn3Utils:Initialize Project`
* copy ``com6RockAvalanche/com6RockAvalancheCfg.ini`` to ``com6RockAvalanche/local_com6RockAvalancheCfg.ini`` and if desired change configuration settings
* if you are on a develop installation, make sure you have an updated compilation, see :ref:`advancedUsage:Update AvaFrame`
* run:
  ::

    python runCom6RockAvalanche.py


Scarp Calculation
-----------------


* first go to ``AvaFrame/avaframe``
* copy ``avaframeCfg.ini`` to ``local_avaframeCfg.ini`` and set your desired avalanche directory name
* create an avalanche directory - for this task you can use :ref:`moduleIn3Utils:Initialize Project`

Input
~~~~~

* All input files are automatically read from the set avalancheDir. No file paths need to be specified
* elevation: DEM (ASCII), which serves as the basis for calculating the scarps. Must be in avalancheDir/Inputs.
* shapefile: A shapefile containing point geometries. These points represent the centers of the ellipsoids or planes.
    The coordinates (x,y) of these points are used. If the plane method is used, the shape file must contain the
    Attributes "zseed", "dip" and "slopeangle" as float values. If the ellipsoid method is used, the shape file must
    contain the Attributes "maxdepth", "semimajor" and "semiminor". The file must be located in avalancheDir/Inputs/
    POINTS and file name must end with “_coordinates”.
* perimeter_shapefile: A shapefile that specifies a boundary area. Must be located in avalancheDir/Inputs/POLYGONS and file name must end with “_perimeter”.

Output
~~~~~~

* elevscarp: Output DGM (ASCII or GeoTIFF), which maps the input DGM minus the calculated scarp. Is saved under
“raster_scarp.asc” in avalancheDir/Inputs/RASTERS.
* hrelease: File path to the output DGM (ASCII or GeoTIFF), which represents the calculated scarp volumes.Is saved
under “raster_rel.asc” in avalancheDir/Inputs/REL.

Config
~~~~~~

Prepare the config file (scarpCfg.ini):

* copy ``com6RockAvalanche/scarpCfg.ini`` to ``com6RockAvalanche/local_scarpCfg.ini`` and if desired change
configuration settings

* Input:
    o set useShapefiles = True
* Settings:
    o method: Here you specify whether the plane or the ellipsoid method should be used

If all the data is provided successfully, start the script by running::

    python runCom6Scarp.py