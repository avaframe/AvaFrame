com3Hybrid: Hybrid modeling
============================

:py:mod:`com3Hybrid` is a computational module that combines the Dense Flow Avalanche (DFA) simulation
model and the Alpha-Beta statistical one and takes advantage
of their strength and try to reduce their weaknesses.
The weakness of the DFA simulation lies in its inputs. It requires, among other a friction parameter :math:`\mu`.
The weakness of the statistical model lies in the path finding (path which will be then used to extract an
avalanche profile and compute a runout angle).
The idea here is to determine the avalanche path in an automated way by running a DFA simulation. The path is
then used to compute the runout angle corresponding to this specific avalanche. The runout is then related to the
friction parameter :math:`\mu` (using a coulomb friction method). This :math:`\mu` is finally used as input parameter to
run a new DFA simulation.


Input
-----

* raster of the DEM (.asc file)
* a release feature (shapefile) in ``Inputs/REL``
* Split point (shapefile). The ultimate goal will be to find the beta point without the need of a split point

Outputs
--------

* The avalanche path specific for the input topography
* the :math:`\mu` value (:math:`\alpha` travel angle) specific for the input topography
* results from the DFA simulation for the input topography
* a combined plot showing the results of the DFA simulation and Alpha-Beta model
* the energy line plot

To run
-------

* go to ``AvaFrame/avaframe``
* copy ``com3Hybrid/com3HybridCfg.ini`` to ``com3Hybrid/local_com3HybridCfg.ini`` and edit (if not, default values are used)
* put your com1DFA and com2AB settings in  ``com3Hybrid/hybridModel_com1DFACfg.ini`` and ``com3Hybrid/hybridModel_com2ABCfg.ini``
  (these files replace the local ``.ini`` files in ``com1DFA`` and ``com2AB``)
* make sure all the required inputs are available in the avalanche directory
* enter the path to the desired dataset in ``local_avaframeCfg.ini``
* run::

      python3 runScripts/runCom3Hybrid.py


Path finding
-------------
The avalanche path is extracted from the DFA simulation.
The first step consists in extracting a path by averaging the trajectory of all particles.
The averaging is weighed by the mass of the particles, it is possible to talk of a mass averaged pass.
The path is then resampled.
The path is extended towards the top of the release. Indeed this input path for the Alpha-Beta model
needs to be consistent with how it was calibrated.



Procedure
-----------

Preprocessing :

* The avalanche path (x,y) is first resampled. Default value for resampling is distance=10m (maximal horizontal distance between two points).
  Note that it does not make much sense to decrease this value to be smaller than the raster grid resolution.
  We then introduce the curvilinear coordinate s which represents the projected horizontal distance along the path.

* The avalanche path is projected on the DEM to generate the profile using a bi-linear interpolation on the DEM to the point of interest.
* The split point (which is not necessarily given on the avalanche path) is projected on the avalanche path.

From this we obtain the (x,y,z) and (s,z) coordinates of the avalanche profile.

AlphaBeta Model:

* Find the :math:`10^\circ` point from (s,z).
* Calculate :math:`\beta`.
* Calculate the :math:`\alpha_j` angles using the adequate standard, small avalanche or custom parameter set.

Postprocessing:

* Plot and save results.

Configuration parameters
---------------------------------

:distance: resampling distance. The given avalanche path is resampled with a step <= 10m (default).

:dsMin: float. Threshold distance [m]. When looking for the beta point make sure at least dsMin meters after the beta point also have an angle bellow 10° (dsMin=30m as default).

:smallAva: boolean (False as default) if True apply :math:`(k_1, k_2, k_3, k_4, SD)` set of small avalanches or False, standard avalanches

:customParam: boolean (False as default). Enables to choose custom :math:`(k_1, k_2, k_3, k_4, SD)``. If True,
  the values from the configuration file are used

:k1: float. Use this value if ``customParam=True``

:k2: float. Use this value if ``customParam=True``

:k3: float. Use this value if ``customParam=True``

:k4: float. Use this value if ``customParam=True``

:SD: float. Use this value if ``customParam=True``


:PlotPath: Plot Avalanche path on raster; default False
:PlotProfile: Plot profile; default False
:SaveProfile: Save profile to file; default True
:WriteRes: Write result to file: default True
