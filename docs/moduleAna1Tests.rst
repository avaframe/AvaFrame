##################################
ana1Tests: Testing
##################################



Dambreak
=========

This function is used to reproduce the exact solution for a dam break problem based on the Savage Hutter model and was described in :cite:`FaMa2012`.
In this test a release mass in suddenly released over an inclined plane.
This function also provides its own plotting routines that can be used to generate plots of the
analytical results and also to compare simulation results to the exact solutions.
An example is given in runDamBreak.py, where the analytical solution is computed and
avalanche simulations are performed and both results are then compared.
The input data for this example can be found in 'data/avaDamBreak'.

The plotResults function generates profile plots for the flow depth and velocity and
an animation of the flow depth evolution.
With the plotComparison function, the simulation results are plotted alongside the
analytical solution for the given time step for both, flow depth and velocity.


Input
-----

damBreakSol:

* avalanche directory
* main avaframe configuration settings
* configuration settings of com1DFA including a section 'DAMBREAK'

plotResults:

* x, h, u - domain extent, flow depth and flow velocity from exact solution
* dtStep - time step (of results)
* cfg - main configuration

plotComparison:

* dataComsol (dataFrame with simulation results)
* hL, xR, hR, uR - initial release thickness, domain extent, flow depth and flow velocity from exact solution
* dtAnalysis - time step
* cfgMain - main configuration

Outputs
-------
damBreakSol:

* hL: initial release thickness
* h: flow depth
* u: flow velocity
* phi: slope angle
* x: horizontal domain

To run
------
An example on how to apply this test to ``com1DFA`` is provided in :py:func:`runScripts/runSimilaritySol`.

* in ``AvaFrame/avaframe`` run::

  runDamBreak.py



Similarity solution
====================

The similarity solution is one of the few cases where an analytic solution can be derived for solving the depth integrated equations.
It is therefore a very useful test case for validating simulation results coming from the dense flow avalanche computation module.
This analytic solution can be derived under very strict conditions and making one major assumption on the shape of the solution.
The full development of the conditions and assumptions as well as the derivation of the solution is presented in details in :cite:`HuSiSaNo1993`.
In this problem, we consider an avalanche governed by a dry friction law (Coulomb friction) flowing down an inclined plane.
The released mass is initially distributed in an ellipse with a parabolic depth shape.
This mass is suddenly released at :math:`t=0` and flows down the inclined plane.

The :py:mod:`ana1Tests.simiSol` module provides functions to compute the analytic solution and some plotting routines
to visualize this solution and to compare it to the output from the DFA computational module.


A workflow example is given in :py:func:`runScripts/runSimilaritySol`, where the analytical solution is computed and
avalanche simulations are performed and both results are then compared.
The input data for this example can be found in ``data/avaSimilaritySol``.

The plotResults function generates profile plots for the flow depth and velocity
in both slope and cross slope directions. The simulation results are plotted alongside the
analytical solution for the given time step.


Input
-----

* avalanche directory (``data/avaSimilaritySol``)
* configuration settings of com1DFA including a section 'SIMISOL' (see ``data/avaSimilaritySol/Inputs/simiSol_com1DFACfg.ini``)


Outputs
-------

solSimi
~~~~~~~~~

Results of the analytic solution of the problem. Dictionary with:

* T: real time array
* f, g and their first derivative: time solution of the transformed variables

getSimiSolParameters
~~~~~~~~~~~~~~~~~~~~~

Function that transforms the f, g transformed variables back to flow depth and
velocity in real time and space.

Plotting routine
~~~~~~~~~~~~~~~~~

Comparing the results from the DFA module to the similarity solution leads to the following plots:


.. list-table::


    * - Along slope flow depth and velocity of the DFA solution and similarity solution

        .. figure:: _static/simiSol_x.png


      - Cross slope flow depth and velocity of the DFA solution and similarity solution

        .. figure:: _static/simiSol_y.png


.. figure:: _static/simiSol_2D.png

    2D flow depth of the DFA solution (full lines) and similarity solution (dashed lines) in meters




Flat plane test
====================

In this test, the behavior of a conical pile of material placed on a flat plane
and ruled by Coulomb friction is studied. The pile starts at rest and depending
on the steepness of the slope and the friction angle will start flowing. The idea behind this test
is to check implementation of the bottom shear stress in the DFA module as well as the gradient
computation (SPH method).

The :py:mod:`ana1Tests.FPtest` module provides functions to compute the analytic solution and some plotting routines
to visualize the pile and to compare it to the output from the DFA computational module.


An workflow example is given in :py:func:`runTestFP`. The input data for this example can be found in ``data/avaFPtest``.


Input
-----

* avalanche directory ('data/avaFPtest')
* configuration settings of com1DFA including a section 'FPSOL' (see ``data/avaFPtest/Inputs/FlatPlane_com1DFACfg.ini``)


Outputs
-------

The :py:func:`plotProfilesFPtest` function generates a profile plot of the flow depth in the radial direction.
The simulation results are plotted alongside the analytical solution for the given time step.


.. list-table::


    * - Pile at initial time step

        .. figure:: _static/flatPlaneTest.png


      - Pile after 19.8s

        .. figure:: _static/flatPlaneTest20s.png
