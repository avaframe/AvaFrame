com1DFA: Module DFA-Kernel
===========================

This module is a simulation tool for dense flow (snow) avalanches (DFA). It is based on the
samos-AT (Snow Avalanche Modeling and  Simulation- Advanced Technologies) developed by the Austrian government
in cooperation with the company AVL List GmbH in Graz.
The calculation of the DFA is based on the depth integrated governing equations and
solved numerically using the smoothed particle hydrodynamic (sph) method.

Dense flow avalanche simulations can be performed for different release area scenarios, with or without
entrainment and/or resistance areas.
There is the option to vary the internal friction parameter or the release snow thickness.


Input
---------

The module requires an avalanche directory, that follows a specified folder structure.
This avalanche directory can be created by running: ``runInitializeProject.py``.
In the directory Inputs, the following files are required:

* digital elevation model as .asc file
* release area as shapefile (in Inputs/REL)

and the following files are optional:

* entrainment area as shapefile (in Inputs/ENT)
* resistance area as shapefile (in Inputs/RES)

The simulation settings area defined in the configuration file ``com1DFACfg.ini``:

* flagOut - print full model output
* flagEnt - include entrainment areas
* flagRes - include resistance areas
* flagVarPar - perform parameter variation
* flagVarEnt - perform parameter variation with entrainment and resistance
* varPar - parameter to be varied
* varParValues - values for parameter variation


Output
---------
The simulation results are saved to: *Outputs/com1DFA* and include:

* raster files of the peak values for pressure, flow depth and speed (*Outputs/com1DFA/peakFiles*)
* reports of all simulations (*Outputs/com1DFA/reports*)
* log files of all simulations
* experiment log that lists all simulations


To run
--------

* create an avalanche directory with required input files - for this task you can use :ref:`moduleIn3Utils:Initialize Project`
* copy ``avaframeCfg.ini`` to ``local_avaframeCfg.ini`` and set your desired avalanche directory name 
* copy ``com1DFACfg.ini`` to ``local_com1DFACfg.ini`` and set your desired settings
* in ``avaframe/`` run:

    python3 runCom1DFA.py


Theory
--------


The governing equations of the dense flow avalanche are derived from the
incompressible mass and momentum balance on a Lagrange control volume
([Zw2000]_ [ZwKlSa2003]_). Assuming the avalanche is much longer and larger
than thick, it is possible to integrate the governing equations over the thickness
of the avalanche and operate some simplifications due to the shape of the avalanche.
This leads, after some calculation steps described in details in Appendix
:ref:`appendixCom1DFA:Governing Equations for the Dense Flow Avalanche` to:

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
See appendix :ref:`appendixCom1DFA:Numerics` for further details.
