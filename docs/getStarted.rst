Getting started
================

Even-though most of the routines in Avaframe are automatized, running them still
requires a minimum amount of python scripting.
Here we provide an example on how to perform simulations with Avaframe,
but first some information on the prerequisites.


Prerequisites
---------------

When running Avaframe, the python module `logging <https://docs.python.org/3/library/logging.config.html>`_ is used.

.. include:: moduleLogging.rst


In order to set the configurations required by all the modules within Avaframe, the python module
`configparser <https://docs.python.org/3/library/configparser.html>`_.


.. include:: moduleConfiguration.rst


How to run - an example
------------------------


The following example should make
it easier for you to find your way in Avaframe.
Before starting, make sure you have installed Avaframe and the required
dependencies as explained in :ref:`installation:Installation`.

Initialize project
~~~~~~~~~~~~~~~~~~
The next step is to create the folder where the input data lies and where the
output results will be saved. To do so, specify the full path to the folder
you want to create in the ``local_avaframeCfg.ini`` (which is a copy of
``avaframeCfg.ini`` that you need to create). Then run ::

  python path/to/runInitializeProject.py

This will create a new directory with the input required by Avaframe structured as described
in :ref:`moduleIn3Utils:Initialize Project`. You are now ready to start.

Input data
~~~~~~~~~~~~

Check the input data required by the different modules you want to use and fill
the ``Inputs/`` folder accordingly. For example for the com1DFA module you nee to
supply the inputs described in :ref:`moduleCom1DFA:Input`. You can also have a
look at the default setting for the module you want to use (for example
``com1DFACfg.ini`` for module com1DFA). You can change those settings in a
``local_`` copy of the ``.ini`` file.

Building your run script
~~~~~~~~~~~~~~~~~~~~~~~~~

In this section, we give an example of the workflow to run Avaframe.
You can then create your own workflow taking the ``runOperational.py`` script as
template. Start by creating a python file (``myAvaframeRunScript.py`` for example).

Imports
~~~~~~~~

The run script starts with importing python modules ::

  # Load python modules
  import os
  import time
  import logging

Then you need to import the avaframe modules you will be using. This list will
grow as you use more avaframe modules. You can refer to the different
computational modules documentation to know which imports are required::

  # Local avaframe imports
  # first import the always useful tools
  from avaframe.in3Utils import cfgUtils
  from avaframe.in3Utils import logUtils
  from avaframe.in3Utils import initializeProject as initProj
  from avaframe.log2Report import generateReport as gR
  # then depending on which computational module you want to use
  # for the Com1DFA module
  from avaframe.com1DFA import com1DFA
  from avaframe.out3SimpPlot import outPlotAllPeak as oP
  # for the Alpha Beta module
  from avaframe.com2AB import com2AB
  from avaframe.out3Plot import outAB

General configuration
~~~~~~~~~~~~~~~~~~~~~

The next step is to get the general settings::

  # Load avalanche directory from general configuration file
  # v0.1 -> one avalanche is fine,
  # v0.2 -> if given dir contains more than one avalancheDir -> run loop
  cfgMain = cfgUtils.getGeneralConfig()
  avalancheDir = cfgMain['MAIN']['avalancheDir']


Logging
~~~~~~~

Then you need to prepare for logging::

  # log file name; leave empty to use default runLog.log
  logName = 'runOperational'
  # Start logging
  log = logUtils.initiateLogger(avalancheDir, logName)
  log.info('MAIN SCRIPT')
  log.info('Current avalanche: %s', avalancheDir)

If you just created the ``avalancheDir`` this one should be clean but if you
already did some calculations you might want to clean it::

  # Clean input directory(ies) of old work and output files
  initProj.cleanSingleAvaDir(avalancheDir, keep=logName)

Computation Modules
~~~~~~~~~~~~~~~~~~~

Now starts the real calls to the computation modules.
You can if necessary time the different routines::

  # Time the routine
  startTime = time.time()

Followed by, at the end::

  # Print time needed
  endTime = time.time()
  log.info('Took %s seconds to calculate.' % (endTime - startTime))

If you want to run the Com1DFA routine, you can do the following::

  # Load input parameters from configuration file
  cfgCom1DFA = cfgUtils.getModuleConfig(com1DFA)
  # Run dense flow
  reportDictList = com1DFA.com1DFAMain(cfgCom1DFA, avalancheDir)
  # Collect results/plots/report  to a single directory
  # Generate plots for all peakFiles
  plotDict = oP.plotAllPeakFields(avalancheDir, cfgCom1DFA, cfgMain['FLAGS'])

If you want to run the Com2AB routine, you would want to do::

  # Load input parameters from configuration file
  cfgAB = cfgUtils.getModuleConfig(com2AB)
  # Run Alpha Beta
  resAB = com2AB.com2ABMain(cfgAB, avalancheDir)
  # Collect results/plots/report  to a single directory
  reportDictList, _, _ = outAB.writeABpostOut(resAB, cfgAB, reportDictList)


Report
~~~~~~

It is also possible to generate a report file (:ref:`moduleLog2Report:Generate Report`)
in addition to the ``.log`` file::

  # Set directory for report
  reportDir = os.path.join(avalancheDir, 'Outputs')
  # write report
  gR.writeReport(reportDir, reportDictList, cfgMain['FLAGS'], plotDict)
