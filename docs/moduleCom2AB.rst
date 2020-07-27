com2AB: Module AlphaBeta 
==========================

Calculates the run-out of an avalanche for a given DEM, avalanche path and split point according to the alphabeta model.

Input
-----

* raster of the DEM
* avalanche path (set of x y coordinates of the path)
* x y coordinates of the split point

Outputs:
--------

* profile plot with alpha, beta and run-out points
* txt file with angle and coordinates and of the different points

To run:
-------

**Deprecated**

* cp local_settings_plot.py
* enter paths to the desired datasets and output files
* python3 AlphaBeta.py


Theory
~~~~~~

Example equation:

.. math::

    A_i(z) = X_i^{wet}(z) - \tau ^{*} \, max \left( T(z) - T_{bulb}, 0 \right)

