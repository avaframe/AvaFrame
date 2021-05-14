Dense Flow Avalanche algorithm and workflow
============================================

Initialization:
-----------------

Read and prepare DEM
~~~~~~~~~~~~~~~~~~~~~~~~
Read and check DEM (remesh if needed).
Prepare DEM for simulation (compute surface normals vector field, cell area)
done in the :py:func:`com1DFAPy.com1DFA.initializeMesh` function:

.. code-block:: python

  # -----------------------
    # Initialize mesh
    demOri, dem = initializeMesh(cfgGen, demOri, methodMeshNormal)

Start of iteration n:
----------------------

Add artificial viscosity
~~~~~~~~~~~~~~~~~~~~~~~~

Compute friction forces
~~~~~~~~~~~~~~~~~~~~~~~~

Compute the bottom shear force according to the friction model chosen and the
resistance force if activated.


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


Take gravity and lateral pressure forces into account
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Take friction into account
~~~~~~~~~~~~~~~~~~~~~~~~~

.. graphviz:: com1DFAAlgorithmGraph.dot
