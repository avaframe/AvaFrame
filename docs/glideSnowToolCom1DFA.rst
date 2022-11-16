Glide snow tool com1DFA
========================

The com1DFA computational module provides the option to add an elastic cohesion force between particles. This can
be used for example to simulate small glide snow avalanches.

Input
-------

The following parameters are required for activating and configuring cohesion:

  * one activates the cohesion by setting ``cohesion`` to 1

  * the maximum strain before breaking of the bond ``cohesionMaxStrain``

  * the Young modulus in N/m² ``cohesiveSurfaceTension``

Theory
--------

Initialization
~~~~~~~~~~~~~~~

During the initialization process, the particles are bonded to each other. This defines the connectivity between
the particles as well as the length at rest of the bonds between the particles.

The construction of the bonds is done by building a triangular mesh on the particles (used as a point cloud).

Cohesion force
~~~~~~~~~~~~~~~

During the computation, an elastic cohesion force is added to the total forces. The force of particle l on particle k
(assuming particle k and l are bonded) reads:

.. math::

  \mathbf{F}_{lk} = -\mathbf{F}_{kl} =  \mathbf{x_\text{old}} - 
