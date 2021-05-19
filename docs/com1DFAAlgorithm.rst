Dense Flow Avalanche algorithm and workflow
============================================

Initialization:
-----------------
At the beginning of the simulation, the avalanche folder and the configuration
are read (:ref:`configuration:Configuration`).
Input data is read and checked according to the chose configuration.
Mesh, particles and fields are sequently initialized.

Initialize Mesh
~~~~~~~~~~~~~~~~~
Read DEM ascii file provided in the Input folder (one and only one should be provided).
If the DEM cell size is different from the :``meshCellSize`` specified in the configuration
from more then ``meshCellSizeThreshold`` [m] the DEM is remeshed (:py:func:`in3Trans.geoTrans.remesh`).
Prepare DEM for simulation, compute surface normals vector field, cell area (:ref:`DFAnumerics:Mesh`). This is done
in the :py:func:`com1DFAPy.com1DFA.initializeMesh` function.

Go back to :ref:`com1DFAAlgorithm:Algorithm graph`

Initialize release, entrainment and resistance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Read and check shape file according to the configuration (check consistency between
what is required by the configuration file and what is available in the ``Inputs`` folder).
Convert shape file features (polygons) to rasters (:py:func:`com1DFAPy.com1DFA.prepareAreas`).
Check consistency of rasters according to the following rules:

  - multiple release features in the release and secondary release shape files
    are allowed but they should not overlap. If they do, simulation terminates with
    an error message.

  - secondary release and entrainment rasters should no overlap between each other or with the
    main release. If they do, the overlapping part is removed. Order of priority is: main
    release, secondary release, entrainment area.

Initialize particles
~~~~~~~~~~~~~~~~~~~~~
Particles are initialize according to the release raster extracted from the release shape file
and the mass per particle determination method (``massPerParticleDeterminationMethod``) specified in the configuration.
The mass per particle determination method can be chosen between:

  - MPPDIR= mass per particle direct. The ``massPerPart`` value is taken from the configration

  - MPPDH= mass per particles through release thickness. The ``massPerPart`` value is computed
    using the release thickness per particle ``deltaTh`` value given in the configration, the area of
    the release mesh cell: :math:`massPerPart = \rho cellArea deltaTh`.

The number of particles placed in each release cell is computed according to the ``massPerPart``
and the area and release thickness of the cell. The number should be an integer meaning that the
float is rounded up or down with a probability corresponding to the decimal part (i.e. 5.7 will
be rounded to 6 with a probability of 0.7 and 5 with a probability of 0.3). This ensures a better
match with the desired ``massPerPart`` value. Particles are then place randomly within the
mesh cell.
Other particles properties velocity, cell number... are also initialized here.
See :py:func:`com1DFAPy.com1DFA.initializeParticles`

Initialize fields
~~~~~~~~~~~~~~~~~
All fields (grid values defined as a raster) ar initialized. Flow velocity, pressure, peak flow velocity and peak pressures
are set to zero. Flow depth and peak flow depth are set according to the initial particle distribution.
See :py:func:`com1DFAPy.com1DFA.initializeFields`


Time scheme and iterations:
------------------------------
The mass and momentum equations described in :ref:`theoryCom1DFA:Governing Equations for the Dense Flow Avalanche` are solved numerically
in time using an operator splitting method. The different forces involved are sequently added to update the velocity.
Position is then updated using a centered Euler scheme.
The time step can either be fixed or dynamically computed using the Courant–Friedrichs–Lewy (CFL) condition.


Compute Forces:
-----------------
This section gives an overview of the different steps to compute the forces acting on the snow particles.
Those forces are separated in several terms: A gravity driving fore (:math:`F_{drive}`), a friction force
(:math:`F_{fric}`), an entrainment force (related to the entrained mass of snow) and an artificial viscous force.
Those forces are computed by the two following functions
:py:func:`com1DFAPy.DFAfunctionsCython.computeForceC` and :py:func:`com1DFAPy.DFAfunctionsCython.computeForceSPHC`

Artificial viscosity
~~~~~~~~~~~~~~~~~~~~~~

In :ref:`theoryCom1DFA:Governing Equations for the Dense Flow Avalanche`, the governing
equations for the DFA were derived and all first order or smaller terms where neglected.
Among those terms is the lateral shear stress. This term leads toward
the homogenization of the velocity field. It means that two neighbor elements
of fluid should have similar velocities. The aim behind adding artificial viscosity is to
take this phenomena into account. The following vicosity force is added:


.. math::
    \begin{aligned}
    \mathbf{F_{viscosity}} = &- \frac{1}{2}\rho C_{Lat}\|\mathbf{du}\|^2 A_{Lat}
    \frac{\mathbf{du}}{\|\mathbf{du}\|}\\
    = & - \frac{1}{2}\rho C_{Lat}\|\mathbf{du}\| A_{Lat} \mathbf{du}
    \end{aligned}

Where the velocity difference reads :math:`\mathbf{du} = \mathbf{u} - \mathbf{\bar{u}}`
(:math:`\mathbf{\bar{u}}` is the mesh velocity interpolated at the particle position).
:math:`C_{Lat}` is a coefficient that rules the viscous force. It would be the
equivalent of :math:`C_{Drag}` in the case of the drag force. The :math:`C_{Lat}`
is a numerical parameter that depends on the mesh size. Its value is set to 100
and should be discussed and further tested.

Adding the viscous force
"""""""""""""""""""""""""

The viscous force is added implicitly:

.. math::
  \begin{aligned}
  \mathbf{F_{viscosity}} = &-\frac{1}{2}\rho C_{Lat}\|\mathbf{du}^{old}\| A_{Lat}
  \mathbf{du}^{new}\\
  = &  -\frac{1}{2}\rho C_{Lat}\|\mathbf{u}^{old} - \mathbf{\bar{u}}^{old}\| A_{Lat}
  (\mathbf{u}^{new} - \mathbf{\bar{u}}^{old})
  \end{aligned}

Updating the velocity is done in two steps. First adding the explcit term related to the
mean mesh velocity and then the implicit term which leads to:

.. math::
  \mathbf{u}^{new} = \frac{\mathbf{u}^{old} - C_{vis}\mathbf{\bar{u}}^{old}}{1 + C_{vis}}

With :math:`C_{vis} = \frac{1}{2}\rho C_{Lat}\|\mathbf{du}^{old}\| A_{Lat}\frac{dt}{m}`


Compute friction forces
~~~~~~~~~~~~~~~~~~~~~~~~
The friction force encompasses all forces that oppose the motion of the particles.
More details about One of those forces is the bottom shear force. The other is an optional resistance force.
Both components are added to the :math:`F_{fric}` force term.

Bottom shear force
"""""""""""""""""""""
This force accounts for the friction between the snow particles and the bottom surface (:ref:`theoryCom1DFA:Bottom friction`).
The expression of the bottom shear stress depends on the friction model chosen but can be written in the
following general forme, :math:`\tau^{(b)}_i = f(\sigma^{(b)},\overline{u},\overline{h},\rho_0,t,\mathbf{x})`.
The friction model and its parameters can be set in the configuration file. More details about the different
friction models are given in :ref:`theoryCom1DFA:Friction model`.
Be aware that the normal stress on the bottom surface :math:`\sigma^{(b)}` is composed of the normal component of the
gravity force and the curvature acceleration term as shown in :eq:`sigmab`. It is possible
to deactivate the curvature acceleration component of the shear stress by setting the
``curvAcceleration`` coefficient to 0 in the configuration file.

Added resistance force
"""""""""""""""""""""""
An additional friction force called resistance can be added. This force aims to model the added
resistance due to the specificity of the terrain on which the avalanche evolves, for example
due to forests. To add a resistance force, one must provide a resistance shape file in the ``Inputs``
folder and switch the ``simType`` to ``res``, ``entres`` or ``available`` to take this resistance area into account.
Then, during the simulation, all particles flowing through this resistance area will undergo an
extra resistance force. More details about how this force is computed and the different parameters chosen
are found in :ref:`Resistance <theoryCom1DFA:Resistance:>`.


Compute body driving force
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This force takes into account the gravity force, which is the driving force of the snow motion.
The expression of this force is rater simple, it represents the tangential (tangent to the surface) part of the gravity force
(the normal part of the force is accounted for in the friction term).


Take entrainment into account
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the entrainment areas, particles can entrain mass through erosion or plowing process.
In both mechanisms, one must account for three things:

    - The change of mass due to the entrainment.

    - The change of momentum. Indeed, the entrained snow was accelerated from rest to the speed of the avalanche.

    - The loss of momentum due to the plowing or erosion phenomena. The entrained mass bounds with the ground
    needs to be broken.

These 3 terms are further detailed in :ref:`Entrainment <theoryCom1DFA:Entrainment:>`.

In the numerics, the mass is updated according to the entrainment model in
:py:func:`com1DFAPy.DFAfunctionsCython.computeEntMassAndForce`. The velocity is updated immediately
after using an implicit formulation.

Compute lateral pressure forces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The lateral pressure forces are related to the gradient of the flow depth. This gradient
is computed using a smoothed particle hydrodynamic method.
This force is added to the :math:`F_{SPH}`.


Update position
----------------

Driving force, lateral pressure force and friction forces are sequently used to update the velocity.
Then the particle position is updated using a centered Euler scheme.
This steps are done in :py:func:`com1DFAPy.DFAfunctionsCython.updatePositionC`.

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
The particles position is updated using the new velocity and a centered Euler scheme.

.. math::
  \mathbf{}^{new} = \mathbf{}^{old} + dt * 0.5 * (\mathbf{u}^{old} + \mathbf{u}^{new})


Correction step:
~~~~~~~~~~~~~~~~

Add secondary release area
----------------------------

Update fields
--------------

Update fields
~~~~~~~~~~~~~

Update particles flow depth
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Algorithm graph
----------------

.. graphviz:: com1DFAAlgorithmGraph.dot
