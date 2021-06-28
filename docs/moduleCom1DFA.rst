com1DFA: DFA-Kernel
===========================

This module is a simulation tool for dense flow (snow) avalanches (DFA). It is based on
samosAT (Snow Avalanche Modeling and  Simulation- Advanced Technologies) developed by the Austrian government
in cooperation with the company AVL List GmbH in Graz.
The calculation of the DFA is based on the depth integrated governing equations and
solved numerically using the smoothed particle hydrodynamics (sph) method.

Dense flow avalanche simulations can be performed for different release area scenarios, with or without
entrainment and/or resistance areas.
The model configuration is controlled via a configuration file.
This configuration file can be modified in order to change any of the default settings and also allows
to perform simulations for varying parameters all at once.


Input
---------

DFA simulations are performed within an avalanche directory that has to be organized following a specified folder structure.

.. Note::  An avalanche directory can be created by running: ``runInitializeProject.py``, which creates the required folder structure:

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
* release area scenario as shapefile (in Inputs/REL); multiple are possible
  -> the release area name should not contain an underscore, if so '_AF' is added

and the following files are optional:

* entrainment area as shapefile (in Inputs/ENT)
* resistance area as shapefile (in Inputs/RES)
* secondary release area as shapefile (in Inputs/SECREL)


Model configuration
--------------------
The model configuration is read from a configuration file: ``com1DFACfg.ini``. In this file,
all model parameters are listed and can be modified. We recommend to create a local copy of the file,
and keep the default configuration in ``com1DFACfg.ini`` untouched.
For this purpose, in ``avaframe/`` run:
  ::

      cp com1DFA/com1DFACfg.ini com1DFA/local_com1DFACfg.ini

and modify the parameter values in there.

Another option is to directly provide the path to a particular configuration file.
The order is as follows, read configuration file from provided path, read local configuration file,
if not present read default configuration file.

It is also possible to perform multiple simulations at once, with varying input parameters.
There are multiple options to vary a parameter:

* replace the default parameter value with desired value
* provide a number of parameter values separated by ``|`` (e.g. ``relTh=1.|2.|3.``)
* provide a number of parameter values using ``start:stop:numberOfSteps`` (e.g. ``relTh=1.:3.:3``)


Output
---------
Using the default configuration, the simulation results are saved to: *Outputs/com1DFA* and include:

* raster files of the peak values for pressure, flow depth and flow velocity (*Outputs/com1DFA/peakFiles*)
* raster files of the peak values for pressure, flow depth and flow velocity for the initial time step (*Outputs/com1DFA/peakFiles/timeSteps*)
* markdown report including figures for all simulations (*Outputs/com1DFA/reports*)
* mass log files of all simulations (*Outputs/com1DFA*)
* configuration files for all simulations (*Outputs/com1DFA/configurationFiles*)


However, in the configuration file, it is possible to change the result parameters and time Steps that shall be exported.
Have a look at the designated Output section in ``com1DFACfg.ini``.


To run
--------

* first go to ``AvaFrame/avaframe``
* copy ``avaframeCfg.ini`` to ``local_avaframeCfg.ini`` and set your desired avalanche directory name
* create an avalanche directory with required input files - for this task you can use :ref:`moduleIn3Utils:Initialize Project`
* copy ``com1DFACfg.ini`` to ``local_com1DFACfg.ini`` and if desired change configuration settings
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
See theory :ref:`theoryCom1DFA:Numerics` for further details.
