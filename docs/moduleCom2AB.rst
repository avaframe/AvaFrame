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

* cp local_com2ABCfg.py
* enter paths to the desired datasets and output files
* in AvaFrame/avaframe/ run python3 runCom2AB.py


Theory:
~~~~~~

.. math::

    \alpha_j = k_1 \beta + k_2 f" + k_3 H_0 + k_4 + j SD

where k_1, k_2, k_3 and k_4 are coefficients and SD the standard deviation
Index j={-1,-2,0,1} and \alpha_j= \alpha + j SD

Numerics:
~~~~~~~~
Pre-processing :
* The avalanche path (x,y) is first re-sampled. Default value for re-sampling is distance=10m (maximal horizontal distance between two points).
  We introduce the curvilinear coordinate s.
* The avalanche path is projected on the DEM to generate the profile using the closest lower left corner of the DEM to the point of interest.
* The split point is projected on the avalanche path.
We end up with the (x,y,z) and (s,z) coordinates of the avalanche profile

AlphaBeta Model:
* Find the 10° point from (s,z).
* Calculate \beta.
* Calculate the \alpha_j angles using the adequate standard ot small avalanche parameter set.

Post-processing:
* Plot and save results.

Optional configuring parameters:
~~~~~~
* distance: re-sampling distance. The given avalanche path is re-sampled with a 10m (default) step.
* smallAva: is True or False (default) depending on if you want to apply the (k_1, k_2, k_3, k_4, SD) set of small avalanches or standard avalanche
* customParam: enables to choose an custom (k_1, k_2, k_3, k_4, SD). customParam = None as default. Otherwise customParam has to be declared as a python dictionary:
customParam = {}
customParam['k1'] = your value
customParam['k2'] = your value
customParam['k3'] = your value
customParam['k4'] = your value
customParam['SD'] = your value
