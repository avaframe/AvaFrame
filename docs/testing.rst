##################################
Testing
##################################


In AvaFrame we provide a continuously growing test suite. This includes tests that target model verification and model validation.
The former are used to test the numerical implementation of the mathematical model, whereas the latter are aiming at testing if
the mathematical model employed is appropriate to describe the physical processes of the desired application, in our case the simulation of dense snow avalanches.

In the following, a brief description of the available tests is given.

.. Note::  This section is currently under development and will be updated regularly.

.. Note::  See this `example pdf <https://avaframe.org/pdf/standardTestsReportExample.pdf>`_ (25mb)
           how we apply these tests during development.

Tests for model verification
=============================

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
Here, we provide an implementation of one of the seven cases presented in :cite:`FaMa2012`.
In Test 2, the initial velocity is zero and a granular mass is suddenly released and flows downslope.
For further details have a look at the module :py:mod:`ana1Tests.damBreak` and the Section :ref:`moduleana1Tests:Dambreak`.



Similarity Solution
--------------------

This test makes use of semi-analytic solution of a sliding granular mass on an inclined plane.
In :cite:`HuSiSaNo1993` the semi-analytic solution is derived for this problem, based on the Savage-Hutter shallow flow
model to describe the granular flow.
Here, we provide an implementation of the similarity solution and compare
it to the output of a numerical simulation.
For further details have a look at the module :py:mod:`ana1Tests.simiSol` and the Section :ref:`moduleana1Tests:Similarity solution`.


Tests for model validation
==========================

In this section, tests that check whether the chosen mathematical model is suitable to represent the problem at hand can be found.
Currently, we provide idealized and real-world topographies that can be used for avalanche simulations. Additionally, there are several functions
that can be used to analyse the simulation results
and compare them to other data such as results from other models or observations. For this, we refer to Sections :ref:`moduleAna3AIMEC:ana3AIMEC: Aimec` and :ref:`moduleOut3Plot:outQuickPlot`.

Idealized test cases
-----------------------

We provide a range of idealized topographies to perform avalanche simulations including different release area scenarios. Some of the topographies also
include entrainment and/or resistance areas. These can be found in ``data/NameOfAvalanche``. The functions to produce these can be found in
module :py:mod:`in3Utils/generateTopo.py`.

Bowl - BL
^^^^^^^^^^

Bowl shaped topopraphy. Used to test e.g. rotational symmetry.

.. figure:: _static/testCaseTopos/DEM_BL_Topo_plot.png


Flat Plane - FP
^^^^^^^^^^^^^^^^

.. figure:: _static/testCaseTopos/myDEM_FP_Topo_plot.png

Inclined Plane - IP
^^^^^^^^^^^^^^^^^^^
IP (Inclined plane)

.. figure:: _static/testCaseTopos/DEM_IP_Topo_plot.png

Parabolic slope - PF
^^^^^^^^^^^^^^^^^^^^
Parabolic slope with flat foreland

.. figure:: _static/testCaseTopos/DEM_PF_Topo_plot.png

Hockeystick - HS
^^^^^^^^^^^^^^^^^

Hockeystick with linear slope and flat foreland and smooth transition

.. figure:: _static/testCaseTopos/DEM_HS_Topo_plot.png

Helix - HX
^^^^^^^^^^^

Helix-shaped topography

.. figure:: _static/testCaseTopos/DEM_HX_Topo_plot.png

Pyramid - PY
^^^^^^^^^^^^^

pyramid-shaped topography, optional with flat foreland

.. figure:: _static/testCaseTopos/DEM_PY_Topo_plot.png

Real-world test cases
----------------------

We provide a range of real-word topographies to perform avalanche simulations including different release area scenarios including
entrainment areas. These can be found in ``data/NameOfAvalanche``.

.. Note::  All images give you a bigger version on click.

Alr
^^^

.. image:: _static/testCaseTopos/avaAlr_plot.png
    :width: 59%
.. image:: _static/testCaseTopos/alrRelEnt.png
    :width: 39%

Gar
^^^

.. image:: _static/testCaseTopos/avaGar_plot.png
    :width: 59%
.. image:: _static/testCaseTopos/garRelEnt.png
    :width: 39%

Hit
^^^

.. image:: _static/testCaseTopos/avaHit_plot.png
    :width: 59%
.. image:: _static/testCaseTopos/hitRelEnt.png
    :width: 39%

Kot
^^^

.. image:: _static/testCaseTopos/avaKot_plot.png
    :width: 59%
.. image:: _static/testCaseTopos/kotRelEnt.png
    :width: 39%

Mal
^^^

.. image:: _static/testCaseTopos/avaMal_plot.png
    :width: 59%
.. image:: _static/testCaseTopos/malRelEnt.png
    :width: 39%

Wog
^^^

.. image:: _static/testCaseTopos/avaWog_plot.png
    :width: 59%
.. image:: _static/testCaseTopos/wogRelEnt.png
    :width: 39%
