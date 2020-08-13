ana3AIMEC: Module Aimec
==========================

Aimec is a post-processing module to analyze results from avalanche simulations.


Input
-----

* raster of the DEM (.asc file)
* avalanche path (as a shape file).
* Results from avalanche simulation

Outputs
--------

* output figures in ``NameOfAvalanche/Outputs/AimecResults/pics/``
* txt file with results in ``NameOfAvalanche/Outputs/AimecResults/``

To run
-------

* copy ``ana3AIMECCfg.py`` to ``local_ana3AIMECCfg.py`` (if not, the demo example is used)
* enter paths to the desired datasets and output files in ``local_ana3AIMECCfg.py`` (if not, the demo example is used)
* in ``AvaFrame/`` run::

      python3 avaframe/runAna3AIMEC.py

Theory
-----------
      .. figure:: _static/aimec_transfo.png
              :width: 90%

              Pressure field on real raster (in blue the chosen path) and on the deskewed raster along given path



Procedure
-----------

Find transformation (from real raster to deskewed raster along given path) :


Projection on results (Speed, Pressure...) on deskewed raster:


Analyze results:

Plot and save results:

Optional configuration parameters
---------------------------------

:domainWidth: width of the domain around the avalanche path in [m]
:pressureLimit: pressure limit value for evaluation of runout in [kPa]
:distance: re-sampling distance. The given avalanche path is re-sampled with a 10m (default) step.
:plotFigure: plot figures; default False
:savePlot: Save figures; default True
:WriteRes: Write result to file: default True
