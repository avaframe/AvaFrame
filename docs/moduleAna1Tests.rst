##################################
ana1Tests: testing
##################################



Dambreak
=========

This function is used to reproduce the exact solution for a dam break problem based on the Savage Hutter model and was described in :cite:`FaMa2012`.
In this test a release mass in suddenly released  over an inclined plane.
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
* configuration settings of com1DFA including a section 'DAMBREAK'

plotResults:

* x, h, u - domain extent, flow depth and flow velocity from exact solution
* dtStep - time step (of results)
* cfg - main configuration

plotComparison:

* dataComsol (dictionary with simulation results)
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

* runDamBreak.py


Similarity solution
====================

This function is used to reproduce the exact similarity solution for a dry friction
avalanche (coulom friction) based on the Savage Hutter model and was described in :cite:`HuSiSaNo1993`.
In this test a release mass (distributed in an elliptic shape) andn suddenly released  over an inclined plane.
This function also provides its own plotting routines that can be used to generate plots of the
analytical results and also to compare simulation results to the exact solutions.
An example is given in runSimilaritySol.py, where the analytical solution is computed and
avalanche simulations are performed and both results are then compared.
The input data for this example can be found in 'data/avaSimilaritySol'.

The plotResults function generates profile plots for the flow depth and velocity
in both slope and cross slope directions. The simulation results are plotted alongside the
analytical solution for the given time step.


Input
-----

* avalanche directory ('data/avaSimilaritySol')
* configuration settings of com1DFA including a section 'SIMISOL'


Outputs
-------
solSimi:

Results of the analytic solution of the problem. Dictionary with:

* T: real time array
* f, g and there first derivative: time solution of the transformed variables

getSimiSolParameters:

Function that transforms the f, g transformed variables back to flow depth and
velocity in real time and space.

To run
------

* runSimilaritySol.py
