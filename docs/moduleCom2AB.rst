com2AB: Module AlphaBeta
==========================

Calculates the run-out of an avalanche for a given DEM, avalanche path and split point according to the alphabeta model.

Input
-----

* raster of the DEM (.asc file)
* set of avalanche paths (as a shape file). There can be multiple paths in the shape file.
* Split point (as a shape file). There can be multiple split points in the shape file.

Outputs:
--------

* profile plot with alpha, beta and run-out points
* txt file with angle and coordinates of the different points

To run:
-------

* copy com2ABCfg.py to local_com2ABCfg.py (if not, the demo example is used)
* enter paths to the desired datasets and output files in local_com2ABCfg.py (if not, the demo example is used)
* in AvaFrame/ run::

      python3 avaframe/runCom2AB.py


Theory:
-------

.. figure:: _static/AlphaBeta_theory.png
        :width: 90%

        Topographical variables for the calculation of :math:`\alpha`


The snow avalanche runout distance is calculated using a statistical model based on data collected for real
avalanches ([Lied_1980]_, [Bakkehøi_1983]_ and [Wagner]_). An equation of the following type is fitted to the data:

.. math::
    \alpha_j = k_1 \beta + k_2 y'' + k_3 H_0 + k_4 + j SD

where :math:`H_0` is the elevation loss of the quadratic fit of the avalanche profile. :math:`y''=\gamma` is the curvature of this same quadratic fit.
:math:`\beta` is the angle of the line between the :math:`10^\circ` point (first point where the avalanche profiles slope is under :math:`10^\circ`)
and the starting point. The coefficients :math:`(k_1, k_2, k_3, k_4)` and the standard deviation :math:`SD` are calculated during the fitting process.
Index :math:`j=\{-1,-2,0,1\}` and :math:`\alpha_j= \alpha + j SD`. These coefficients obviously depend on the initial set of chosen data. :math:`\alpha_0= \alpha`
is the angle between the stopping and the starting point of the avalanche. :math:`\alpha_j= \alpha + j SD` takes into account the uncertainty of the process.

Numerics:
-------

Pre-processing :

* The avalanche path (x,y) is first re-sampled. Default value for re-sampling is distance=10m (maximal horizontal distance between two points). Note that it does not make much sense to decrease this value lower than the raster grid resolution. We then introduce the curvilinear coordinate s.
* The avalanche path is projected on the DEM to generate the profile using a bi-linear interpolation on the DEM to the point of interest.
* The split point (which is not necessarily given on the avalanche path) is projected on the avalanche path.
We end up with the (x,y,z) and (s,z) coordinates of the avalanche profile.

AlphaBeta Model:

* Find the :math:`10^\circ` point from (s,z).
* Calculate :math:`\beta`.
* Calculate the :math:`\alpha_j` angles using the adequate standard, small avalanche or custom parameter set.

Post-processing:

* Plot and save results.

Optional configuration parameters:
-------
* distance: re-sampling distance. The given avalanche path is re-sampled with a 10m (default) step.

::

    distance = 10

* smallAva: is True or False (default) depending on if you want to apply the :math:`(k_1, k_2, k_3, k_4, SD)` set of small avalanches or standard avalanches

::

    smallAva = False

* customParam: enables to choose custom :math:`(k_1, k_2, k_3, k_4, SD)``. customParam = None as default. Otherwise customParam has to be declared as a python dictionary:

::

    customParam = {} # = None by default
    customParam['k1'] = 'your value'
    customParam['k2'] = 'your value'
    customParam['k3'] = 'your value'
    customParam['k4'] = 'your value'
    customParam['SD'] = 'your value'

* Plot save write result flag :

::

  flags = {}
  flags['PlotPath'] = False
  flags['PlotProfile'] = False
  flags['SaveProfile'] = True
  flags['WriteRes'] = True


References
----------

.. [Lied_1980] Lied, K. und Bakkehøi, S. (1980)
    Empirical calculations of snow-avalanche run-out distance based on topographic parameters. In: Journal of Glaciology, Band 26: 165–177.

.. [Bakkehøi_1983] Bakkehøi, S., Domaas, U. und Lied, K. (1983).
    Calculation of snow avalanche runout distance. In: Annals of Glaciology, Band 4: 24–29.

.. [Wagner] Peter Wagne (2016).
    Masterarbet : Kalibrierung des α-β-Modells für das Ermitteln der Auslauflänge von kleinen und mittleren Lawinen
