com2AB: Module AlphaBeta
==========================

Calculates the run-out of an avalanche for a given DEM, avalanche path and split point according to the alphabeta model.

Input
-----

* raster of the DEM (.asc file)
* set if avalanche paths (as a shape file). There can be multiple paths in the shape file.
* Split point (as a shape file). There can be multiple split points in the shape file.

Outputs:
--------

* profile plot with alpha, beta and run-out points
* txt file with angle and coordinates and of the different points

To run:
-------

* cp com2ABCfg.py to local_com2ABCfg.py (if not, the demo example is used)
* enter paths to the desired datasets and output files in local_com2ABCfg.py (if not, the demo example is used)
* in AvaFrame/ run python3 avaframe/runCom2AB.py


Theory:
-------

.. math::
    \alpha_j = k_1 \beta + k_2 f" + k_3 H_0 + k_4 + j SD

where :math:`(k_1, k_2, k_3, k_4)` are coefficients and :math:`SD` the standard deviation
Index :math:`j=\{-1,-2,0,1\}` and :math:`\alpha_j= \alpha + j SD`

Numerics:
-------
Pre-processing :

* The avalanche path (x,y) is first re-sampled. Default value for re-sampling is distance=10m (maximal horizontal distance between two points). We introduce the curvilinear coordinate s.
* The avalanche path is projected on the DEM to generate the profile using the closest lower left corner of the DEM to the point of interest.
* The split point is projected on the avalanche path.
We end up with the (x,y,z) and (s,z) coordinates of the avalanche profile

AlphaBeta Model:

* Find the 10° point from (s,z).
* Calculate :math:`\beta`.
* Calculate the :math:`\alpha_j` angles using the adequate standard of small avalanche parameter set.

Post-processing:

* Plot and save results.

Optional configuring parameters:
-------
* distance: re-sampling distance. The given avalanche path is re-sampled with a 10m (default) step.
* smallAva: is True or False (default) depending on if you want to apply the :math:`(k_1, k_2, k_3, k_4, SD)` set of small avalanches or standard avalanche
* customParam: enables to choose custom :math:`(k_1, k_2, k_3, k_4, SD)``. customParam = None as default. Otherwise customParam has to be declared as a python dictionary:

customParam = {}

customParam['k1'] = your value

customParam['k2'] = your value

customParam['k3'] = your value

customParam['k4'] = your value

customParam['SD'] = your value
