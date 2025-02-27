.. _in4Region:

in4Region: Regional Input Management
====================================
:py:mod:`in4Region` is a module for organizing large amounts of avalanche input data into multiple avalanche directories, based 
on the AvaFrame data structure, thus enabling efficient and controlled input management of larger datasets for further processing.
Importantly, the module also provides a simple, automatic method for clipping the input DEM around each release area group, 
which then forms the basis for the organization of other input. Ultimately, in combination with :ref:`com7Regional <com7Regional>`, 
the module provides an efficient processing workflow for regional scale avalanche simulations with :ref:`com1DFA <com1DFA>`.

.. note::
    This module is still under development and may contain bugs or incomplete features.

Input
-----
The module is currently compatible with the following input file types:

**Required:**
    * Digital elevation model (as .asc or .tif file)
    * Release areas (shapefile in REL directory, with optional attributes ``group`` and ``scenario``)

**Optional:**
    * Entrainment areas (shapefile in ENT directory)
    * Resistance areas (shapefile in RES directory)

The module is intended to work with only one of each of these files. The first file that is found is the one that is used.

Where the expected input directory structure is::

  avalancheDir/
  └── Inputs/
      ├── REL/
      │   └── *.shp         # release areas    
      ├── ENT/              # entrainment areas
      │   └── *.shp
      ├── RES/              # resistance areas
      │   └── *.shp
      └── *.asc or *.tif    # digital elevation model (DEM)

Group and scenario creation
---------------------------
Input data organization is based on two key concepts:

1. **Groups**: Collections of avalanche release areas (single polygon features) that are located in the same spatial domain and may be 
wanted to be simulated together (e.g. divided into avalanche paths)
2. **Scenarios**: Collections of release area features within each group, that would be simulated together in :ref:`com1DFA <com1DFA>`.

These may be defined through two new attributes in the input release area shapefile attribute table:

* ``group``: expected format: string (e.g. *"group1"*)

* ``scenario``: expected format: comma-separated list of strings, without spaces (e.g. *"small,large"* or *"10y,30y,100y"*)

For example:

.. list-table::
    :header-rows: 1
    :widths: 40 40 40

    * - name
      - group
      - scenario
    * - rel1
      - avaPath1
      - 1,all
    * - rel2
      - avaPath1
      - 2,all
    * - rel3
      - avaPath2
      - 1

In this example, four release scenarios (shapefiles) would be created, in two separate avalanche directories::

  avaPath1/
  └── REL/
      ├── avaPath1_1.shp     - containing [1] feature: rel1
      ├── avaPath1_2.shp     - containing [1] feature: rel2
      └── avaPath1_all.shp   - containing [2] features: rel1, rel2

  avaPath2/
  └── REL/
      └── avaPath2_1.shp     - containing [1] feature: rel3

In the case that scenarios are defined for only some release features within a group, the rest will be grouped together as a single 'NULL' scenario.

In the case that no attributes or values for ``group`` or ``scenario`` are provided, the procedure will create groups with single scenario shapefiles for each release area feature, 
which may be wanted to simulate each release feature separately without much effort. Keep in mind however that this approach may result in a large amount of duplicated DEM data,
in the case that release areas are located in close proximity to each other.

Output
------
Running ``runScripts/runSplitInputs.py`` with valid input data will result in the following output in ``<avalancheDir>/Outputs/in4Region/splitInputs``:

1. Individual avalanche directories for each group containing:

  - Clipped DEM file
  - Scenario-specific release area shapefiles
  - Optional: Clipped entrainment and resistance areas

2. Two visual reports (see :numref:`fig-splitInputs-basic` and :numref:`fig-splitInputs-optional`)

3. Scenario report in txt format (see example below)

.. list-table::
   :widths: 50 50

   * - .. _fig-splitInputs-basic:
       .. figure:: /_static/splitInputs_visualReport_basic.png
          :width: 100%
          :alt: Basic visual report

          Example of basic inputs report displaying resulting groups and their extent
     - .. _fig-splitInputs-optional:
       .. figure:: /_static/splitInputs_visualReport_optional.png
          :width: 100%
          :alt: Optional inputs report

          Example of optional inputs report displaying RES and ENT areas for each group

Example Scenario Report:

.. code-block:: text

    SCENARIO REPORT
    ==============
    Generated: 2025-02-04 10:58:56

    Group: group1
    ------------

    Scenario: sce1
    No. of release areas: 2
    - rel1
    - rel2

    No. of entrainment areas: 1
    No. of resistance areas: 3

    Group: group2
    ------------
    ...

Configuration
-------------
Settings are controlled through ``splitInputsCfg.ini``, in which the ``bufferSize`` for the group extent is defined (which is used for DEM, RES, and ENT clipping into 
smaller chunks). By default, this value is set to 2500 m. For each group, a bounding box is created from the maximum x-y extent of all release features in the group. 
The value for ``bufferSize`` is then added to each direction (+x, -x, +y, -y). This buffer may be adjusted according to the expected maximum runout length of your avalanches - 
a larger value will ensure that no simulation will exit its domain, while a larger value will result in smaller output file sizes.

Procedure
---------
The ``splitInputsMain`` function, which is called in ``runScripts/runSplitInputs.py``, performs the following steps:

1. Create central avalanche directory list
2. Set up avalanche directories
3. Split and write release areas to each directory
4. Clip and write DEM to each directory
5. Clip and write optional input to each directory (currently includes RES and ENT)
6. Divide release areas into scenarios
7. Write reports

To Run
------
1. Prepare inputs in your ``<avalancheDir>/Inputs``
2. Configure settings in ``splitInputsCfg.ini`` (or local version ``local_splitInputsCfg.ini``)
3. Set path to avalanche directory in ``avaframeCfg.ini`` (or local version ``local_avaframeCfg.ini``)
4. Execute from AvaFrame/avaframe directory:

.. code-block:: bash

    python3 runScripts/runSplitInputs.py
