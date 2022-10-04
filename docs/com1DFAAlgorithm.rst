com1DFA Algorithm and workflow
==============================

Algorithm graph
----------------

The following graph describes the Dense Flow Avalanche simulation workflow
(the different nodes are clickable and link to the detailed documentation)

.. graphviz:: com1DFAAlgorithmGraph.dot


Initialization:
-----------------
At the beginning of the simulation, the avalanche folder and the configuration
are read (:ref:`configuration:Configuration`).
Input data is fetched according to the chosen configuration.
Mesh, particles and fields are subsequently initialized.

Initialize Mesh
~~~~~~~~~~~~~~~~~

Read DEM ascii file provided in the Input folder (only one DEM ascii file allowed).
If the DEM cell size is different from the :``meshCellSize`` specified in the configuration
from more than ``meshCellSizeThreshold`` [m] the DEM is remeshed (:py:func:`in3Trans.geoTrans.remesh`).

Prepare DEM for simulation, compute surface normals vector field, cell area (:ref:`DFAnumerics:Mesh`).

This is done in the :py:func:`com1DFA.com1DFA.initializeMesh` function.

Go back to :ref:`com1DFAAlgorithm:Algorithm graph`

Initialize release, entrainment and resistance areas
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Read and check shapefiles according to the configuration (check consistency between
what is required by the configuration file and what is available in the ``Inputs`` folder).
Convert shapefile features (polygons) to rasters (:py:func:`com1DFA.com1DFA.prepareArea`).
Check consistency of rasters according to the following rules:

  - multiple release features in the release and secondary release shapefiles
    are allowed but they should not overlap. If they do, simulation terminates with
    an error message.

  - secondary release and entrainment rasters should no overlap between each other or with the
    main release. If they do, the overlapping part is removed. Order of priority is: main
    release, secondary release, entrainment area.

Go back to :ref:`com1DFAAlgorithm:Algorithm graph`

Initialize particles
~~~~~~~~~~~~~~~~~~~~~
Particles are initialized according to the release raster extracted from the release shapefile
and the mass per particle determination method (``massPerParticleDeterminationMethod``) specified in the configuration.
The mass per particle determination method can be chosen between:

  - MPPDIR= mass per particle direct. The ``massPerPart`` value is taken from the configuration and is the same for all
  cells.

  - MPPDH= mass per particles through release thickness. The ``massPerPart`` value is computed
    using the release thickness per particle ``deltaTh`` value given in the configuration and the area of
    the release mesh cell: :math:`\mbox{massPerPart} = \rho\times \mbox{cellArea} \times\mbox{deltaTh}`.

    - MPPKR= mass per particles through number of particles per kernel radius. There is no ``massPerPart`` since it can
    vary from one cell to another depending on the release thickness of the cells. The aim of this method is to ensure a
    constant density of particles within the snow domain (``nPPK`` particles per kernel radius is the target value).
    This is related to the SPH method used for computing the flow thickness gradient. It requires a
    sufficient number of particles to properly approximate the flow thickness gradient. It makes the most sense to combine
    the MPPKR particle initialization method with the `splitOption` 1. In this combination, the particles will be merged
    or split to keep a constant density of particles per kernel radius (:ref:`DFAnumerics:Splitting and merging`).

.. Note::  If MPPDIR is used, consider adapting the mass per particle value when changing the mesh cell size from the default.
           This is important because, when using MPPDIR, the total number of particles is independent of the cell size. Hence,
           reducing the cell size results in less particles per cell, whereas when using MPPDH,
           the number of particles per cell is fixed (considering the respective release thickness and deltaTh value).
           Reducing the cell size will increase the total number of particles but not the number of
           particles per cell. Finally, using the MPPKR method, the number of particles per cell is independent from
           both cell size and release thickness (``nPPK`` particles per kernel radius is the target value).

The number of particles placed in each release cell is computed according to the ``massPerPart`` or ``nPPK`` depending
on the ``massPerParticleDeterminationMethod`` chosen and the area and/or release thickness of the cell.
The number should be an integer meaning that the float is rounded up or down with a probability corresponding to the
decimal part (i.e. 5.7 will be rounded to 6 with a probability of 0.7 and 5 with a probability of 0.3).
This ensures a better match with the desired ``massPerPart`` value. Particles are then placed randomly within the
mesh cell.
Other particles properties velocity, cell number... are also initialized here.
See :py:func:`com1DFA.com1DFA.initializeParticles`.

Go back to :ref:`com1DFAAlgorithm:Algorithm graph`

Initialize fields
~~~~~~~~~~~~~~~~~
All fields (mesh values defined as a raster) are initialized. Flow velocity, pressure, peak flow velocity and peak pressures
are set to zero. Flow thickness and peak flow thickness are set according to the initial particle distribution.
See :py:func:`com1DFA.com1DFA.initializeFields`

Go back to :ref:`com1DFAAlgorithm:Algorithm graph`


Time scheme and iterations:
------------------------------
The mass and momentum equations described in :ref:`theoryCom1DFA:Governing Equations for the Dense Flow Avalanche` are solved numerically
in time using an operator splitting method. The different forces involved are sequentially added to update the velocity
(see :ref:`DFAnumerics:Adding forces`).
Position is then updated using a centered Euler scheme.
The time step can either be fixed or dynamically computed using the Courant–Friedrichs–Lewy (CFL) condition
(in the second case one must set ``cflTimeStepping`` to ``True`` and set the desired CFL coefficient).

Go back to :ref:`com1DFAAlgorithm:Algorithm graph`


Compute Forces:
-----------------
This section gives an overview of the different steps to compute the forces acting on the snow particles.
Those forces are separated in several terms: A gravity driving fore (:math:`F_{drive}`), a friction force
(:math:`F_{fric}`), an entrainment force (related to the entrained mass of snow) and an artificial viscous force.
Those forces are computed by the two following functions
:py:func:`com1DFA.DFAfunctionsCython.computeForceC` and :py:func:`com1DFA.DFAfunctionsCython.computeForceSPHC`.

Go back to :ref:`com1DFAAlgorithm:Algorithm graph`

Artificial viscosity
~~~~~~~~~~~~~~~~~~~~~~

This viscous friction force is artificially added to the numerical computation.
The aim of this force is to stabilize the simulation and prevent neighbor particles
to have too significantly different velocities. Physically, this force also makes sense and corresponds
to some second order forces that were neglected (lateral shear stress) as explained in
:ref:`DFAnumerics:Artificial viscosity`.
This force is controlled by the ``subgridMixingFactor`` in the configuration file.
Setting this parameter to 0 deactivates the artificial viscosity term.
The default value (set to 100) does not have any physical foundation yet. Future work
will help defining this parameter in a more physical way. Remember that the artificial viscosity is dependent on the grid cell size.

The velocity is updated immediately after using an explicit/implicit formulation.

Go back to :ref:`com1DFAAlgorithm:Algorithm graph`



Compute friction forces
~~~~~~~~~~~~~~~~~~~~~~~~
The friction force encompasses all forces that oppose the motion of the particles.
One of those forces is the bottom shear force. The other is an optional resistance force.
Both components are added to the :math:`F_{fric}` force term.

Bottom shear force
"""""""""""""""""""""
This force accounts for the friction between the snow particles and the bottom surface.
The expression of the bottom shear stress depends on the friction model chosen but can be written in the
following general form, :math:`\tau^{(b)}_i = f(\sigma^{(b)},\overline{u},\overline{h},\rho_0,t,\mathbf{x})`.
The friction model is set by the ``frictModel`` value and the corresponding parameters can be set in the configuration file.
More details about the different friction models are given in :ref:`theoryCom1DFA:Friction model`.
Be aware that the normal stress on the bottom surface :math:`\sigma^{(b)}` is composed of the normal component of the
gravity force and the curvature acceleration term as shown in :eq:`sigmab`. It is possible
to deactivate the curvature acceleration component of the shear stress by setting the
``curvAcceleration`` coefficient to 0 in the configuration file.


Added resistance force
"""""""""""""""""""""""
An additional friction force called resistance can be added. This force aims to model the added
resistance due to the specificity of the terrain on which the avalanche evolves, for example
due to forests. To add a resistance force, one must provide a resistance shapefile in the ``Inputs/RES``
folder and switch the ``simType`` to ``res``, ``entres`` or ``available`` to take this resistance area into account.
Then, during the simulation, all particles flowing through this resistance area will undergo an
extra resistance force. More details about how this force is computed and the different parameters chosen
are found in :ref:`Resistance <theoryCom1DFA:Resistance:>`.


Go back to :ref:`com1DFAAlgorithm:Algorithm graph`


Compute driving force
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This force takes into account the gravity force, which is the driving force of the snow motion.
The expression of this force is rather simple, it represents the tangential (tangent to the surface) part of the gravity force
(the normal part of the force is accounted for in the friction term).


Go back to :ref:`com1DFAAlgorithm:Algorithm graph`



Take entrainment into account
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Snow entrainment can be added to the simulation. One must provide an entrainment shapefile
in ``Inputs/ENT`` and set the ``simType`` to ``ent``, ``entres`` or ``available``
(see :ref:`com1DFAAlgorithm:Initialize release, entrainment and resistance areas`).
In the entrainment areas defined by the entrainment shapefile, particles can entrain mass through erosion or plowing.
In both mechanisms, one must account for three things:

    - change of mass due to the entrainment

    - change of momentum - entrained snow was accelerated from rest to the speed of the avalanche

    - loss of momentum due to the plowing or erosion processes -entrained mass bounds with the ground needs to be broken

These three terms are further detailed in :ref:`Entrainment <theoryCom1DFA:Entrainment:>`. The parameters
used to compute these processes can be set in the configuration file.

In the numerics, the mass is updated according to the entrainment model in
:py:func:`com1DFA.DFAfunctionsCython.computeEntMassAndForce`. The velocity is updated immediately
after using an implicit formulation.


Go back to :ref:`com1DFAAlgorithm:Algorithm graph`


Compute lateral pressure forces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The lateral pressure forces (:math:`F_{SPH}`) are related to the gradient of the flow thickness (:ref:`DFAnumerics:Forces discretization`). This gradient
is computed using a smoothed particle hydrodynamic method (:ref:`DFAnumerics:SPH gradient`).


Go back to :ref:`com1DFAAlgorithm:Algorithm graph`

Update position
----------------

Driving force, lateral pressure force and friction forces are subsequently used to update the velocity.
Then the particle position is updated using a centered Euler scheme.
These steps are done in :py:func:`com1DFA.DFAfunctionsCython.updatePositionC`.


Take gravity and lateral pressure forces into account
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:math:`F_{drive}` and :math:`F_{SPH}` are summed up and taken into account to update the velocity.
This is done via an explicit method.

Take friction into account
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:math:`F_{fric}` is taken into account to update the velocity.
This is done via an implicit method.

Update particle position
~~~~~~~~~~~~~~~~~~~~~~~~~
The particles position is updated using the new velocity and a centered Euler scheme:

.. math::
  \mathbf{x}^{new} = \mathbf{x}^{old} + dt 0.5 (\mathbf{u}^{old} + \mathbf{u}^{new})


Correction step:
~~~~~~~~~~~~~~~~
The particles z coordinate it readjusted so that the particles lie on the surface of the slope.
There are two reasons why the particles might not lie on the surface anymore after updating their position
according to the computed velocities:

  - 1) because of the inaccuracy related to the time and space discretization.
    This can lead to a particle position being slightly above or under the surface.
    We want to correct this inaccuracy and therefore reproject the particle on the surface
    using its x and y coordinates.

  - 2) because of the curvature of the slope and the particle velocity, particles can become
    detached from the ground in - in this case, the particle is located above the
    surface. In the current state, the com1DFA kernel does not allow this.
    If a particle becomes detached, the particle is also reprojected onto the surface
    using its x and y coordinates.

Similarly, the particles velocity is corrected in order to ensure that it lies in the tangent
plane to the surface (the velocity vector magnitude is preserved, only the direction is changed).

The way the particles position is reprojected onto the surface does not allow both the
velocity magnitude and the particle displacement to match perfectly. This is amplified
by highly curved topographies or abrupt changes in slope.


Go back to :ref:`com1DFAAlgorithm:Algorithm graph`

Add secondary release area
----------------------------
If a secondary release area is provided, the flow thickness
field from the previous time step is used to release a potential secondary release area.
To do so, the flow thickness field is compared to the secondary release area rasters. If
they overlap, the secondary release area is triggered and the secondary release particles
are initialized and added to the flowing particles.


Go back to :ref:`com1DFAAlgorithm:Algorithm graph`

Update fields
--------------

This steps are done in :py:func:`com1DFA.DFAfunctionsCython.updateFieldsC`.

Update fields
~~~~~~~~~~~~~
The mesh values are updated with the particles properties using
:ref:`particles to mesh interpolation <DFAnumerics:Particles to mesh>` methods.
This is used to compute flow thickness, flow velocity and pressure fields from the particle properties.

Update particles flow thickness
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The mesh flow thickness is finally used to update the particle flow thickness value
using :ref:`mesh to particle interpolation <DFAnumerics:Mesh to particle>` methods.


Go back to :ref:`com1DFAAlgorithm:Algorithm graph`
