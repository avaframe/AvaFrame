.. _com7Regional:

com7Regional: Regional Modeling
===============================
:py:mod:`com7Regional` is an experimental module targeted at the application of the :ref:`com1DFA <com1DFA>` kernel on a regional scale. 
In essence, it is a wrapper module that allows for the concurrent execution of com1DFA within multiple avalanche directories, with the aim of reducing 
the overall computational load through parallelization. The module also provides functions for managing and aggregating outputs after a regional simulation, such as com1DFA peak files and merged rasters.
The log output gives information about the processing status of each individual avalanche directory, as well as a summary of how many simulations were processed successfully.

.. note::
    Experience with :ref:`com1DFA <com1DFA>` is recommended before using this module.

Input
-----
A directory structure containing pre-configured avalanche directories (containing an ``Inputs`` folder) is required. For input preparation, we recommend using 
:ref:`in4Region <in4Region>`, which splits merged input data into standard :ref:`com1DFA <com1DFA>` inputs across multiple avalanche directories.

Example of a valid directory structure::

    regionalDir/
    ├── avalanche1/
    │   └── Inputs/
    │       ├── REL/*.shp
    │       └── *.asc or *.tif
    ├── avalanche2/
    │   └── Inputs/
    │       ├── REL/*.shp
    │       └── *.asc or *.tif      
    └── ...

.. note::
    Note the use of "regionalDir" to distinguish between the main input directory from the lower level avalanche directories. 
    However, the regionalDir is still set in the ``avalancheDir`` parameter within ``avaframeCfg.ini``, or within the command call, like in other AvaFrame modules.

Output
------
Outputs are organized in two levels:

**1. Merged rasters** and **2. Individual outputs (per avalanche directory)**

Merged rasters
^^^^^^^^^^^^^^
Configure in ``com7RegionalCfg.ini`` (or local):

.. code-block:: ini

    [GENERAL]
    mergeOutput = True
    mergeTypes = pfv # Available options: [ppr|pfv|pft|pta|FT|FV|P|FM|Vx|Vy|Vz|TA]
    mergeMethod = max  # Available options: [max|min|mean|sum|count]

Produces merged rasters of all peakFile results found within the avalanche directories, for each ``mergeTypes`` and ``mergeMethod``configured, in
``<regionalDir>/Outputs/com7Regional/mergedRasters/``.

Individual outputs
^^^^^^^^^^^^^^^^^^
After running com7 with the given module (e.g. :ref:`com1DFA <com1DFA>`), the standard output is located within each of the avalanche directories within e.g.
``<regionalDir>/<avalancheDir>/Outputs/com1DFA``. Additionally, com7 provides the option of aggregating all output peakFiles and tSteps results into a single directory 
for easier management, either through copying or moving the files after an executed run.

Configure in ``com7RegionalCfg.ini`` (or local):

.. code-block:: ini

    [GENERAL]
    copyPeakFiles = True
    moveInsteadOfCopy = False

Creates::

    Outputs/com7Regional/
    └── allPeakFiles/
        └── allTimeSteps/

Configuration
-------------
Three configuration files are used (locals are preferred if they exist):

1. Main configuration (``avaframeCfg.ini``)
    - To set nCPUs for handling the amount of avalanche directories processed in parallel
    - Handle plot and other output generation
    
2. com7Regional configuration (``com7RegionalCfg.ini``)
    - Manages output aggregation and merged raster creation
    - Overrides com1DFA parameters if specified

3. com1DFA configuration (``com1DFACfg.ini``)
    - Standard simulation parameters

Processing
----------
Parallelization is handled through the concurrent.futures library, specifically the 
`ProcessPoolExecutor class <https://docs.python.org/3/library/concurrent.futures.html#concurrent.futures.ProcessPoolExecutor>`_. 
In essence, tasks are executed concurrently within each input avalanche directory, based on the number of currently available CPUs. 
The maximum number of CPUs is set by the ``nCPU`` parameter in ``avaframeCfg.ini``. By default, to avoid nested parallelization, each 
avalanche directory is assigned a single CPU , essentially meaning that any variations (e.g. through different 
scenarios, parameter variations, etc.) within each avalanche directory are handled sequentially. As a consequence, if the number of 
variations is high, and the number of avaDirs to process is lower than ``nCPU``, it may be more efficient to run simulations with the 
standard ``runCom1DFA.py`` instead, to utilize its parallel processing of variations. Alternatively, advanced users may want to adjust the nCPU for 
variations in ``com7Regional.py``.

To Run
------
1. Prepare input directories, we recommend using :ref:`in4Region <in4Region>`
2. Configure settings in ``com7RegionalCfg.ini`` (or local version ``local_com7RegionalCfg.ini``)
3. Set path to regional directory in ``avaframeCfg.ini`` (or local version ``local_avaframeCfg.ini``)
4. Execute from AvaFrame/avaframe directory:

.. code-block:: bash

    python3 runScripts/runCom7Regional.py
