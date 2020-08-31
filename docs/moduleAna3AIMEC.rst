ana3AIMEC: Module Aimec
==========================

Aimec is a post-processing module to analyze and compare results from avalanche simulations.
It enables the comparison of different simulations (with different input parameters variations for example)
of the same avalanche (meaning using the same DEM and going down the same avalanche path) in a standardized way.


Inputs
-------

* raster of the DEM (.asc file)
* avalanche path in LINES (as a shape file named ``path_aimec.shp``).
* a splitPoint in POINTS (as a shape file named ``splitPoint.shp``).
* Results from avalanche simulation (when using results from com1DFA,
  the helper function ``mainDfa2Aimec`` in ``dfa2Aimec.py`` fetches and prepares the input for Aimec)

Outputs
--------

* output figures in ``NameOfAvalanche/Outputs/ana3AIMEC/com1DFA/pics/``
* txt file with results in ``NameOfAvalanche/Outputs/ana3AIMEC/com1DFA/``
  (a detailed list of the results is described in :ref:`analyze-results`)

To run
-------

* copy ``ana3AIMECCfg.py`` to ``local_ana3AIMECCfg.py`` (if not, the standard settings are used)
* enter path to the desired ``NameOfAvalanche/`` folder in ``AvaFrame/avaframe/avaframeCfg.ini``
* in ``AvaFrame/avaframe/`` run::

      python3 runAna3AIMEC.py

.. _Theory:

Theory
-----------

The simulation results (two dimensional fields of e.g. peak pressure or flow depth) are processed in a way
that it is possible to compare characteristic values such as for example, run-out, maximum peak pressure or
flow depth, maximum peak velocity or speed for different simulations.


AIMEC (Automated Indicator based Model Evaluation and Comparison, [Fischer2013]_) was developed
to analyze and compare avalanche simulations. The computational module presented here is inspired from the original AIMEC code.
The simulations are analyzed and compared by projecting the results along a chosen poly-line (same line for all the simulations
that are compared) called avalanche path.
The raster data, initially located on a regular and uniform grid (with coordinates x and y) is projected on a regular non uniform grid
(grid points are not uniformly spaced) that follows the avalanche path (with curvilinear coordinates (s,l)).
This grid can then be "straightened" or "deskewed" in order to plot it in the (s,l) coordinates system.
The following figure illustrates the process.

.. list-table::



    * - .. _fig-aimec-comp-real:

        .. figure:: _static/aimec_comparison_real_topo.png

            In the real coordinate system (x,y)

      - .. _fig-aimec-comp-new:

        .. figure:: _static/aimec_comparison_new_topo.png

           In the new coordinate system (s,l)


.. _mean-max-path:

Mean and max values along path
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
All two dimensional field results (for example peak pressure or flow depth) can be projected into the curvilinear system using
the previously described method. The maximum and average values of those fields are computed in each cross-section (l direction).
For example the maximum :math:`A_{cross}^{max}(s)` and average :math:`\bar{A}_{cross}(s)` of the two
dimensional distribution :math:`A(s,l)` is:

.. math::
    A_{cross}^{max}(s) = \max_{\forall l \in [-\frac{w}{2},\frac{w}{2}]} A(s,l) \quad\mbox{and}\quad
    \bar{A}_{cross}(s) = \frac{1}{w}\int_{-\frac{w}{2}}^{\frac{w}{2}} A(s,l)dl

Run-out point
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The run-out point corresponding to a given pressure threshold :math:`P_{lim}>0kPa` is the first point :math:`s=s_{runout}`
where the maximum peak pressure falls below the pressure limit (:math:`P_{cross}^{max}(s)<P_{Lim}`). This :math:`s=s_{runout}` is related
to a :math:`(x_{runout},y_{runout})` in the original coordinate system. It is very important to note that the position of this
point depends on the chosen pressure limit value. It would also be possible to use :math:`\bar{P}_{cross}(s)<P_{Lim}` instead of
:math:`P_{cross}^{max}(s)<P_{Lim}`.

Run-out length
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This length depends on what is considered the beginning of the avalanche :math:`s=s_{start}`. It can be related to the release area,
to the transition point (first point where the slope angle is below :math:`30^{\circ}`) or to the run-out area point
(first point where the slope angle is below :math:`10^{\circ}`). The run-out length is then defined as :math:`L=s_{runout}-s_{start}`.

Mean and max indicators
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
From the maximum values along path of the distribution :math:`A(s,l)` calculated in :ref:`mean-max-path`, it is possible to calculate
the global maximum (MMA) and average maximum (AMA) values of the two dimensional distribution :math:`A(s,l)`:

.. math::
    MMA = \max_{\forall s \in [s_{start},s_{runout}]} A_{cross}^{max}(s) \quad\mbox{and}\quad
    AMA = \frac{1}{s_{runout}-s_{start}}\int_{s_{start}}^{s_{runout}} A_{cross}^{max}(s)ds


Area indicators
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
When comparing the run-out area (corresponding to a given pressure threshold) of two simulations,
it is possible to distinguish four different zones. For example, if the first simulation (sim1) is taken as reference and if True corresponds
to the assertion that the avalanche covered this zone and False there was no avalanche in this zone, those four zones are:

    * TP (true positive) zone: green zone on :numref:`fig-aimec-comp-new` , sim1 = True  sim2 = True
    * FP (false positive) zone: blue zone on :numref:`fig-aimec-comp-new` , sim1 = False  sim2 = True
    * FN (false negative) zone: red zone on :numref:`fig-aimec-comp-new` , sim1 = True  sim2 = False
    * TN (true negative) zone: gray zone on :numref:`fig-aimec-comp-new` , sim1 = False  sim2 = False

The two simulations are identical (in the run-out zone) when the area of both FP and FN is zero. In order to provide a normalized
number describing the difference between two simulations, the area of the different zones is normalized by the area of the reference
simulation :math:`A_{ref} = A_{TP} + A_{FP}`. This leads to the 4 area indicators:

    * :math:`\alpha_{TP} = A_{TP}/A_{ref}`, which is 1 if sim2 covers at least the reference
    * :math:`\alpha_{FP} = A_{FP}/A_{ref}`, which is a positive value if sim2 covers an area outside of the reference
    * :math:`\alpha_{FN} = A_{FN}/A_{ref}`, which is a positive value if the reference covers an area outside of sim2
    * :math:`\alpha_{TN} = A_{TN}/A_{ref}`

Identical simulations (in the run-out zone) lead to :math:`\alpha_{TP} = 1` , :math:`\alpha_{FP} = 0` and :math:`\alpha_{FN} = 0`

Mass indicators
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
From the analysis of the release mass (:math:`m_r` at the beginning, i.e :math:`t = t_{ini}`), total mass
(:math:`m_t` at the end, i.e :math:`t = t_{end}`) and entrained mass (:math:`m_e` at the end, i.e :math:`t = t_{end}`)
it is possible to calculate the growth index :math:`GI` and growth gradient :math:`GG` of the avalanche:

.. math::
    GI = \frac{m_t}{m_r} = \frac{m_r + m_e}{m_r} \quad\mbox{and}\quad GG = \frac{m_r + m_e}{t_{end}-t_{ini}}

Procedure
-----------

This section describes how the theory is implemented in the ``ana3AIMEC`` module.

Make Domain transformation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
First, the transformation from (x,y) coordinate system (where the original rasters lie in) to (s,l) coordinate system is applied
given a new domain width. A new grid corresponding to the new domain (following the avalanche path) is built.
The transformation information are stored in a ``rasterTransfo`` dictionary:

:xllc: x coordinate of the lower left cell of the (x,y) domain
:yllc: y coordinate of the lower left cell of the (x,y) domain
:cellsize: original raster cell size
:domainWidth: desired width for the new domain
:gridx: x coordinate of the new raster points (2D numpy array of size (n,m))
:gridy: y coordinate of the new raster points (2D numpy array of size (n,m))
:s: new s coordinates (1D numpy array of size n)
:l: new l coordinates (1D numpy array  of size m)
:x: x coordinate of the centerline (s,l=0) of the new raster (1D numpy arrayof size n)
:y: y coordinate of the centerline (s,l=0) of the new raster (1D numpy arrayof size m)
:rasterArea: area of the cells of the new raster grid (2D numpy array of size (n,m))
:indSplit: index of the projected split point on the avalanche path
:runoutAngle: run-out angle value (in degres)
:indRunoutPoint: index of the run-out point (first point under the given runoutAngle)

Assign data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The simulation results (for example peak pressure or flow depth) are projected on the new grid using the
transformation information. The projected results are stored in the ``newRasters`` dictionary.

.. _analyze-results:

Analyze results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Calculates the different indicators described in the :ref:`Theory` section for a given pressure threshold.
Returns a ``resAnalysis`` dictionary with the analysis results.

:runout: (x,y) coordinates of the run-out as well as the run-out length based on P_cross_max and the pressure Threshold
:runoutMean: (x,y) coordinates of the run-out as well as the run-out length based on P_cross_mean and the pressure Threshold
:AMPP: average maximum peak pressure
:MMPP: maximum maximum peak pressure
:AMD: average maximum flow depth
:MMD: maximum maximum flow depth
:elevRel: z coordinate of the release area (first point with max Peak pressure over pressure Threshold)
:deltaH: DeltaZ between the release point and run-out point
:relMass: release Mass
:entMass: entrained Mass
:growthIndex: growth Index
:growthGrad: growth Gradient
:pressureLimit: pressure Threshold
:pCrossAll: :math:`P_{cross}^{max}(s)` for each simulation

.. _plot-save-results:

Plot and save results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Plots and saves the desired figures. Writes results in ``resAnalysis`` to a text file.
By default, Aimec saves five plots plus as many plots as numerical simulations to
compare to the reference. The first five ones are :
  * "domTransfo" shows the real domain on the left and new domain on the right
  * "referenceFields" shows the peak pressure, flow depth and speed in the new domain
  * "slComparisom" shows the difference the statistics in terms of the peak pressure along the avalanche of all the simulations
  * "ROC" shows the normalized area difference between reference and other simulations.
  * "rmmp" shows the relative difference in maximum peak pressure between reference and other simulation function of run-out length
The last plots "_i_compToRef" where "i" gives the number of the simulation plots the 2D difference with the reference.

Configuration parameters
---------------------------------

:domainWidth: width of the domain around the avalanche path in [m]
:pressureLimit: pressure limit value for evaluation of runout in [kPa]
:distance: re-sampling distance. The given avalanche path is re-sampled with a 10m (default) step.
:plotFigure: plot figures; default False
:savePlot: Save figures; default True
:WriteRes: Write result to file: default True


References
----------

.. [Fischer2013] Fischer, Jan-Thomas. (2013).
    A novel approach to evaluate and compare computational snow avalanche simulation. Natural Hazards and Earth System Sciences. 13. 1655-. 10.5194/nhess-13-1655-2013.
