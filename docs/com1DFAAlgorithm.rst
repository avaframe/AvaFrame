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

Initialize fields
~~~~~~~~~~~~~~~~~


Time scheme and iterations:
------------------------------



Compute Forces:
-----------------

Artificial viscosity
~~~~~~~~~~~~~~~~~~~~~~~~

Compute friction forces
~~~~~~~~~~~~~~~~~~~~~~~~

Bottom shear force
"""""""""""""""""""""

Compute the bottom shear force according to the friction model chosen and the
resistance force if activated.


Added resistance force
"""""""""""""""""""""""

Compute body driving force
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Compute the gravity force component.


Take entrainment into account
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Update mass according to the entrainment model.
Update velocity (momentum conservation and dissipation)


Compute lateral pressure forces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Compute SPH gradients that lead to the lateral pressure forces.


Update position
----------------

Take gravity and lateral pressure forces into account
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Take friction into account
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Update particle position
~~~~~~~~~~~~~~~~~~~~~~~~~

Add secondary release area
----------------------------

Update fields
--------------

Update fields
~~~~~~~~~~~~~

Update particles flow depth
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. graphviz:: com1DFAAlgorithmGraph.dot
