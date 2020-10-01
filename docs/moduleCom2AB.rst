com2AB: Module AlphaBeta
==========================

Calculates the run-out of an avalanche for a given DEM, avalanche path and split point according to the alphabeta model
suited for avalanches in the Austrian Alps.

Input
-----

* raster of the DEM (.asc file)
* set of avalanche paths (as a shape file). There can be multiple paths in the shape file.
* Split point (as a shape file). There can be multiple split points in the shape file.

Outputs
--------

* profile plot with alpha, beta and run-out points
* txt file with angle and coordinates of the different points

To run
-------

* copy ``com2ABCfg.py`` to ``local_com2ABCfg.py`` (if not, the demo example is used)
* enter paths to the desired datasets and output files in ``local_com2ABCfg.py`` (if not, the demo example is used)
* in ``AvaFrame/`` run::

      python3 avaframe/runCom2AB.py


Theory
------

.. figure:: _static/AlphaBeta_theory.png
        :width: 90%

        Topographical variables for the calculation of :math:`\alpha`


The snow avalanche run-out distance is calculated using a statistical model based on data collected for real
avalanches (:cite:`LiBa1980,BaDoLi1983,Wa2016`). An equation of the following type is fitted to the data:

.. math::
    \alpha_j = k_1 \beta + k_2 y'' + k_3 H_0 + k_4 + j SD

where :math:`H_0` is the elevation loss of the quadratic fit of the avalanche profile. :math:`y''=\gamma` is the curvature of this same quadratic fit.
:math:`\beta` is the angle of the line between the :math:`10^\circ` point (first point where the avalanche profiles slope is under :math:`10^\circ`)
and the starting point. The coefficients :math:`(k_1, k_2, k_3, k_4)` and the standard deviation :math:`SD` are calculated during the fitting process.
Index :math:`j=\{-1,-2,0,1\}` and :math:`\alpha_j= \alpha + j SD`. These coefficients obviously depend on the initial set of chosen data. :math:`\alpha_0= \alpha`
is the angle between the stopping and the starting point of the avalanche. :math:`\alpha_j= \alpha + j SD` takes into account the variability of the process.
The values of the SD used are based on normal distribution. It is important to note that a bigger run-out angle leads to a shorter run-out distance. This means
that :math:`\alpha_{-1}= \alpha - SD` leads to a longer run-out. In other words, the probability of the run-out being shorter
than :math:`x_{\alpha_{-1}}` corresponding to :math:`\alpha_{-1}` is approximately 83%.


In this module, the coefficients :math:`(k_1, k_2, k_3, k_4)` and the standard deviation :math:`SD` are already known, they are simply used in the :math:`\alpha`
equation to calculate the run-out on a new profile.

Procedure
-----------

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

Optional configuration parameters
---------------------------------

:distance: re-sampling distance. The given avalanche path is re-sampled with a 10m (default) step.

:smallAva: is True or False (default) depending on if you want to apply the :math:`(k_1, k_2, k_3, k_4, SD)` set of small avalanches or standard avalanches

:customParam: enables to choose custom :math:`(k_1, k_2, k_3, k_4, SD)``. customParam = None as default. Otherwise customParam has to be declared as a python dictionary:

::

    customParam = {} # = None by default
    customParam['k1'] = 'your value'
    customParam['k2'] = 'your value'
    customParam['k3'] = 'your value'
    customParam['k4'] = 'your value'
    customParam['SD'] = 'your value'

:PlotPath: Plot Avalanche path on raster; default False
:PlotProfile: Plot profile; default False
:SaveProfile: Save profile to file; default True
:WriteRes: Write result to file: default True
