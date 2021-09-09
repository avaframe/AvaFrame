Frequently Asked Questions
==========================

Here, we collect questions that we found helpful for users who want to develop their own simulation workflow
or develop new code.


What about changing the spatial resolution of the simulations preformed with com1DFA?
-----------------------------------------------------------
The default setup of :py:mod:`com1DFA` has been tested (add refs.) to produce reasonable results for large
dense flow snow avalanches (size range).
For this reason, if a new simulation is started with a DEM input of cell size other than 5 meters,
com1DFA first remeshes the provided input DEM to 5 meters and then, the computations are performed using the remeshed DEM.
There is also the option to use a different mesh cell size, however keep in mind that there is no tested setup available
for this case.
Still, we want to list here a few things that are useful to take into account when changing the mesh cell size.


When changing the *meshCellSize* in the configuration file of :py:mod:`com1DFA` (``com1DFACfg.ini`` or rather your local copy of it),
consider also changing the *sphKernelRadius* to match the new mesh cell size. Regarding the computation of the
gradients (:ref:`DFAnumerics:SPH gradient`), which are required to compute the lateral forces,
besides the *sphKernelRadius* also the number of particles per cell has
to adjusted - if the *sphKernelRadius* is smaller (due to a higher resolution of the input DEM and hence mesh cell size),
also the number of particles has to be increased in order to ensure a reasonable estimate of the gradient using the SPH
method. If you are using the default setup for computing the mass per particle (*massPerParticleDeterminationMethod*
-see :ref:`com1DFAAlgorithm:Initialize particles`),
which is *MPPDH* (mass per particle through release thickness), the number of particles per cell is independent of
the mesh cell size and hence, does not have to be adjusted when changing the mesh cell size.
However, if *MPPDIR* (mass per particle set directly) is chosen, the
number of particles per cell depends on the release snow mass within a cell, which is computed using the cell area and
the release thickness - hence, the number of particles per cell is a function of the mesh cell size. In this case,
if the mesh cell size and the *sphKernelRadius* are decreased, the *massPerParticle* value should be decreased
accordingly, in order to ensure that a reasonable amount of particles are used to compute the gradients.
