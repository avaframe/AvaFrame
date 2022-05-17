com1DFA: DFA-Kernel
===========================

:py:mod:`com1DFA` is a module for dense flow (snow) avalanche computations (DFA) .
It is a python and cython implementation of the DFA C++ implementation samosAT
(Snow Avalanche Modeling and  Simulation- Advanced Technologies) developed by the Austrian government
in cooperation with the company AVL List GmbH in Graz (see :ref:`moduleCom1DFAOrig:com1DFAOrig: Original DFA-Kernel`).
Calculations are based on the thickness integrated governing equations and
solved numerically using the smoothed particle hydrodynamics (sph) method. Please note 
the use of *thickness averaged/integrated* instead of *depth averaged/integrated* for clarity and consistency. 

Dense flow avalanche simulations can be performed for different release area scenarios, with or without
entrainment and/or resistance areas, and is controlled via a configuration file.
The configration can be modified in order to change any of the default settings and also allows
to perform simulations for varying parameters all at once.


Input
---------

DFA simulations are performed within an avalanche directory, organized with the
folder structure described below.

.. Note::  An avalanche directory can be created by running: :py:mod:`runInitializeProject.py`, which creates the required folder structure:

  ::

    NameOfAvalanche/
      Inputs/
        REL/      - release area scenario
        RES/      - resistance areas
        ENT/      - entrainment areas
        POINTS/   - split points
        LINES/    - avalanche paths
        SECREL/   - secondary release areas
      Outputs/
      Work/


In the directory ``Inputs``, the following files are required:

* digital elevation model as .asc file
  -> use `ESRI grid format <https://desktop.arcgis.com/en/arcmap/10.3/manage-data/raster-and-images/esri-ascii-raster-format.htm>`_
* release area scenario as shapefile (in Inputs/REL); multiple features are possible
  -> the release area name should not contain an underscore, if so '_AF' is added

and the following files are optional:

* entrainment area as shapefile (in Inputs/ENT)
* resistance area as shapefile (in Inputs/RES)
* secondary release area as shapefile (in Inputs/SECREL)

Release-, entrainment thickness settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Release, entrainment and secondary release thickness can be specified in two different ways:

1. Via **shape file**:
  - add an attribute called `thickness` for each feature
  - important: ALL features have to have a single thickness value, which can differ between features
  - for backwards compatibility, the attribute 'd0' also works, but we suggest to use `thickness` in new projects
  - set the flag `THICKNESSFromShp` (i.e. relThFromShp, entThFromShp, secondaryRelthFromShp) to True in the configuration file (default is True)
  - a parameter variation can be added with the `THICKNESSPercentVariation` parameter in the configuration file in the
    form of ``+-percentage$numberOfSteps``. Provided a `+` a positive variation will be performed,
    if `-` is given, only a negative variation is performed. If no sign is given: both directions will be used.
    Additionally, a variation can be added with the `THICKNESSRangeVariation` parameter in the configuration file in the
    form of ``+-range$numberOfSteps``. Provided a `+` a positive variation will be performed,
    if `-` is given, only a negative variation is performed. If no sign is given: both directions will be used.

2. Via **configuration file (ini)**:
  - set the flag 'THICKNESSFromShp' to False
  - provide your desired thickness value in the respective THICKNESS parameter (i.e. relTh, entTh or secondaryRelth)
  - in addition to the `THICKNESSPercentVariation` and `THICKNESSRangeVariation` options (see option 1) and the standard variation options in
    :ref:`configuration:Configuration`, you can also directly set e.g. `relTh = 1.$50$2`, ``referenceValue$+-percentage$numberOfSteps``,
    resulting in a variation of relTh from 0.5 to 1.5m in two steps.

Only available for release thickness:

3. Via **release thickness file**:
  - set the flag 'relThFromShp' to False
  - set the flag 'relThFromFile' to True
  - save a raster file with info on release thickness as .asc file in ``Inputs/RELTH``
    the number of rows and columns must match the DEM raster with desired meshCellSize



DEM input data
^^^^^^^^^^^^^^^^
Regarding the DEM data: if the DEM in ``Inputs`` is not of cell size 5 meters, it is remeshed to a
cell size of 5 meters. However, it is also possible to specify a desired cell size in the
configuration file (parameter `meshCellSize`). In this case, also consider reading :ref:`FAQ:Can the spatial resolution of simulations performed with com1DFA (dense flow) be changed?`.
If the cell size of the DEM in ``Inputs`` is equal to the desired mesh cell size, the DEM is used without modification. If the cell sizes do not match, several options are available:

    - cleanDEMremeshed = True, directory ``Inputs/DEMremeshed`` is cleaned, and the DEM in Inputs/
      is remeshed to the desired cell size - this is the default setting

    - cleanDEMremeshed = False and a DEM including the name of the DEM in Inputs/ and the desired cell size is found
      in Inputs/DEMremeshed - this DEM is used without modification

    - cleanDEMremeshed = False and no matching DEM is found in Inputs/DEMremeshed - the DEM in Inputs/ is remeshed
      to the desired cell size

If the DEM in Inputs/ is remeshed, it is then saved to ``Inputs/DEMremeshed` and available for subsequent
simulations.


Model configuration
--------------------
The model configuration is read from a configuration file: ``com1DFA/com1DFACfg.ini``. In this file,
all model parameters are listed and can be modified. We recommend to create a local copy
and keep the default configuration in ``com1DFA/com1DFACfg.ini`` untouched.
For this purpose, in ``AvaFrame/avaframe/`` run:

  ::

    cp com1DFA/com1DFACfg.ini com1DFA/local_com1DFACfg.ini

and modify the parameter values in there. For more information see :ref:`configuration:Configuration`.

It is also possible to perform multiple simulations at once, with varying input parameters.


Output
---------
Using the default configuration, the simulation results are saved to: *Outputs/com1DFA* and include:

* raster files of the peak values for pressure, flow thickness and flow velocity (*Outputs/com1DFA/peakFiles*)
* raster files of the peak values for pressure, flow thickness and flow velocity for the initial time step (*Outputs/com1DFA/peakFiles/timeSteps*)
* markdown report including figures for all simulations (*Outputs/com1DFA/reports*)
* mass log files of all simulations (*Outputs/com1DFA*)
* configuration files for all simulations (*Outputs/com1DFA/configurationFiles*)


However, in the configuration file, it is possible to change the result parameters and time Steps that shall be exported.
Have a look at the designated subsection Output in ``com1DFA/com1DFACfg.ini``.


To run
--------

* first go to ``AvaFrame/avaframe``
* copy ``avaframeCfg.ini`` to ``local_avaframeCfg.ini`` and set your desired avalanche directory name
* create an avalanche directory with required input files - for this task you can use :ref:`moduleIn3Utils:Initialize Project`
* copy ``com1DFA/com1DFACfg.ini`` to ``com1DFA/local_com1DFACfg.ini`` and if desired change configuration settings
* if you are on a develop installation, make sure you have an updated compilation, see
  :ref:`installation:Setup AvaFrame`
* run:
  ::

    python3 runCom1DFA.py


Theory
--------


The governing equations of the dense flow avalanche are derived from the
incompressible mass and momentum balance on a Lagrange control volume
([Zw2000]_ [ZwKlSa2003]_). Assuming the avalanche is much longer and larger
than thick, it is possible to integrate the governing equations over the thickness
of the avalanche and operate some simplifications due to the shape of the avalanche.
This leads, after some calculation steps described in details in Theory
:ref:`theoryCom1DFA:Governing Equations for the Dense Flow Avalanche` to:

.. math::
    \begin{aligned}
    &\frac{\mathrm{d}V(t)}{\mathrm{d}t} = \frac{\mathrm{d}(A_b\overline{h})}{\mathrm{d}t}
    = \frac{\rho_{\text{ent}}}{\rho_0}\,w_f\,h_{\text{ent}}\,\left\Vert \overline{\mathbf{u}}\right\Vert\\
    &\frac{\,\mathrm{d}\overline{u}_i}{\,\mathrm{d}t} =
    g_i + \frac{K_{(i)}}{\overline{\rho}\,A\,\overline{h}}\,\oint\limits_{\partial{A}}\left(\frac{\overline{h}\,\sigma^{(b)}}{2}\right)n_i\,\mathrm{d}l
    -\delta_{i1}\frac{\tau^{(b)}}{\overline{\rho}\,\overline{h}} - C_{\text{res}}\,\overline{\mathbf{u}}^2\,\frac{\overline{u_i}}{\|\overline{\mathbf{u}}\|}
    -\frac{\overline{u_i}}{A\,\overline{h}}\frac{\,\mathrm{d}(A\,\overline{h})}{\,\mathrm{d}t} + \frac{F_i^{\text{ent}}}{\overline{\rho}\,A\,\overline{h}}\\
    &\overline{\sigma}^{(b)}_{33} = \rho\,\left(g_3-\overline{u_1}^2\,\frac{\partial^2{b}}{\partial{x_1^2}}\right)\,\overline{h}
    \end{aligned}


Numerics
---------

Those equations are solved numerically using a **SPH** method (:cite:`LiLi2010,Sa2007`).
**SPH**  is a mesh free method where the basic idea is to divide the avalanche into
small mass particles. The particles interact with each other according to the
equation of motion described in :ref:`moduleCom1DFA:Theory` and the chosen kernel function.
This kernel function describes the domain of influence of a particle (through the smoothing length parameter).
See theory :ref:`theoryCom1DFA:com1DFA DFA-Kernel theory` for further details.
