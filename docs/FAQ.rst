Frequently Asked Questions
==========================

.. Note::
   If you want to develop your own simulation workflow or develop new code we also suggest to have a look at the code
   examples provided in the `runScripts` directory.


Can the spatial resolution of simulations performed with com1DFA (dense flow) be changed?
-----------------------------------------------------------

The default setup of :py:mod:`com1DFA` is calibrated for medium to large dense flow snow avalanches in regards to hazard
mapping. To provide calibrated results, new simulations with DEM cell sizes other than 5 meters are being
remeshed by com1DFA to **5 meters** for computation.

It is possible to change the cell size, but keep in mind there are no tested/calibrated setups available.
We want to list a few things that we strongly suggest to take into account when changing the mesh cell size:

- Start by changing the *meshCellSize* in the configuration file of :py:mod:`com1DFA` (``com1DFACfg.ini`` or rather your
  local copy of it), and consider changing the *sphKernelRadius* to match the new mesh cell size.

- For the :ref:`DFAnumerics:SPH gradient`, required to compute the lateral forces, the number of particles
  per cell has to be adjusted. If the *sphKernelRadius* (see first point) is decreased, the number of particles should
  increase to ensure a reasonable estimate of the gradient. There are three options for setting the number of particles
  which is computing from the mass per particle
  (*massPerParticleDeterminationMethod*; see :ref:`com1DFAAlgorithm:Initialize particles`):

  - With the default setup *MPPDH* (mass per particle through release thickness), the number of particles per cell is
    independent of the mesh cell size. So no adjustment is necessary.

  - If *MPPKR* (mass per particles through number of particles per kernel radius) is chosen, the size of particles is
    adjusted to get a constant number of particles within an *sphKernelRadius*. This ensures a reasonable amount of
    particles for the gradient computation.

  - If *MPPDIR* (mass per particle set directly) is chosen, the number of particles per cell depends on the release snow
    mass within a cell. The number of particles is computed using the cell area and the release thickness. To ensure a
    reasonable amount of particles for decreased mesh size and *sphKernelRadius*, the *massPerParticle* value should be
    decreased accordingly. This ensures a reasonable amount of particles for the gradient computation.

Regarding DEM cell sizes: 
- if a cell size XX different than the default 5 meters is chosen 
- and the provided DEM in Inputs is of different mesh cell size than XX
- the folder `DEMremeshed` will be checked for a matching DEM.
- if no matching DEM is found -> remeshing takes place

Using the default settings, the subdirectory 'DEMremeshed' is cleaned when starting simulations, to ensure remeshing. This can be changed in the configuration file (further details are provided in :ref:`moduleCom1DFA:DEM input data`).
