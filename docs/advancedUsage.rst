Advanced Usage (Script)
=======================

After the installation, the next sections describe how to get started. The sections configuration and logging
describe the general methods we use, this is helpful to understand how you can change model parameters and similar. 

.. include:: developinstall.rst

--------------------------

.. include:: developinstallwin.rst

--------------------------

.. include:: gettingstarted.rst

---------------------------

.. include:: configuration.rst

-------------------------

.. include:: logging.rst

-------------------------

Example runscripts
------------------

In :py:mod:`runScripts` we provide ready-to-use scripts for different applications of the modules provided within
AvaFrame.


Derive input data
^^^^^^^^^^^^^^^^^^^
- :py:mod:`runScripts.runComputeDist`

Create a new project
^^^^^^^^^^^^^^^^^^^^^^
- :py:mod:`runScripts.runInitializeProject`


Generate idealized/generic topography data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- :py:mod:`runScripts.runGenerateTopo`
- :py:mod:`runScripts.runGenProjTopoRelease`


Postprocessing
^^^^^^^^^^^^^^^

- :py:mod:`runScripts.runAna3AIMEC`
- :py:mod:`runScripts.runAna3AIMECCompMods`
- :py:mod:`runScripts.runAna3AIMECIndi`
- :py:mod:`runScripts.runStatsExample`
- :py:mod:`runScripts.runProbAna`


Visualisation
^^^^^^^^^^^^^^^
- :py:mod:`runScripts.runQuickPlotSimple`
- :py:mod:`runScripts.runQuickPlotOne`
- :py:mod:`runScripts.runPlotTopo`
- :py:mod:`runScripts.runExportToCsv`


Testing
^^^^^^^^
- :py:mod:`runScripts.runDamBreak`
- :py:mod:`runScripts.runSimilaritySol`
- :py:mod:`runScripts.runTestFP`
- :py:mod:`runScripts.runStandardTestsCom1DFAOrig`
- :py:mod:`runScripts.runComparisonModules`
- :py:mod:`runScripts.runFetchBench`
- :py:mod:`runScripts.runWriteDesDict`

------------------------------

Update AvaFrame
------------------------------

To update go to your avaframe repository [YOURDIR]/Avaframe,  pull the latest changes via::

  git pull

and repeat the compilation step from above::

  python setup.py build_ext --inplace

If there are updates on the requirements inside ``setup.py``, it might be necessary to run::

  pip install -e .

again to get the additional requirements installed.