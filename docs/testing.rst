##################################
Testing
##################################


In AvaFrame we attempt to provide a continuously growing test suite. This includes test that strive towards model verification and model validation.
The former are used to test the numerical implementation of the mathematical model, whereas the latter are aiming at testing if
the mathematical model employed is appropriate to describe the physical processes of the desired application, in our case the simulation of dense snow avalanches.

In the following, a brief description of the available tests is given.

.. Note::  This section is currently under development and will be updated regularly.

Tests for model verification
=========================

This section includes testing of bits and pieces of the model code up to tests that check whether the full workflow is numerically robust.
However, at the 'smallest scale', where individual functions are tested, we also refer to the Section about *pytests*: :ref:`develop:How to test code`.


Rotational Symmetry
--------------------
In this test, a pyramid-shaped topography is used to test whether different grid alignments in respect to the topography and hence flow direction cause problems for the model implementation.
For this test, we simply provide two pyramid shaped topographies, where one is rotated around the z-axis to change the orientation of the pyramid's facets with respect to the mesh.

Dambreak problem
------------------

This test is based on a Riemann problem, where the initial condition is described by two different states that are separated by a discontinuity.
In :cite:`FaMa2012` exact solutions are derived for this problem, based on the Savage-Hutter model to describe the granular flow.
Here, we provide an implementation of two of the seven cases presented in :cite:`FaMa2012`.
In Test 2, the initial velocity is zero and a granular mass is suddenly released and flows downslope, whereas in Test 3 the release thickness is constant over the domain, however there is a
discontinuity in velocity. The upper part is has an initial velocity of the same magnitude, however in reverse direction, i.e. introducing a separation of the release mass.
For further details have a look at the module ``ana1Test`` and the Section :ref:`moduleana1Tests:Dambreak`.



Similarity Solution
--------------------

This test makes use of an exact solution of a sliding granular mass on an inclined plane.
In :cite:`HuSiSaNo1993` the analytic solutions is derived for this problem, based on the Savage-Hutter shallow flow
model to describe the granular flow.
Here, we provide an implementation of the analytic similarity solution and compare
it to the output of a numerical simulation.
For further details have a look at the module ``ana1Test`` and the Section :ref:`moduleana1Tests:Similarity solution`.


Tests for model validation
========================

In this section, tests that check whether the chosen mathematical model is suitable to represent the problem at hand. Currently, we provide idealized and real-world Topographies
that can be used for avalanche simulations. Additionally, there are several functions that can be used to analyse the simulation results
and compare them to other data such as results from other models or observations. For this, we refer to Sections :ref:`moduleAna3AIMEC:ana3AIMEC: Aimec` and :ref:`moduleOut3Plot:quickPlot`.

Idealized test cases
-----------------------

We provide a range of idealized topographies to perform avalanche simulations including different release area scenarios. Some of the topographies also
include entrainment and/or resistance areas. These can be found in ``data/NameOfAvalanche``.

Real-world test cases
----------------------

We provide a range of real-word topographies to perform avalanche simulations including different release area scenarios including
entrainment areas. These can be found in ``data/NameOfAvalanche``.
