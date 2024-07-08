com6RockAvalanche: Rock Avalanche
=================================

.. Warning:: This is highly experimental and not tested!

The com6RockAvalanche computational module provides an override setting for com1DFA targeting the simulation of rock
avalanches.

Input
-------

The standard inputs required to perform a simulation run using :py:mod:`com1DFA` 
can be found here: :ref:`moduleCom1DFA:Input`.
However there is one main difference: com6RockAvalanche NEEDS a release thickness raster file. This file has to have
the exact same dimensions as the topography file.
There is a run script to perform a rock avalanche com1DFA run: :py:mod:`runCom6RockAvalanche.py`,
and the configuration settings can be found in ``com6RockAvalanche/com6RockAvalancheCfg.ini``.

To run
------

* first go to ``AvaFrame/avaframe``
* copy ``avaframeCfg.ini`` to ``local_avaframeCfg.ini`` and set your desired avalanche directory name
* create an avalanche directory with required input files - for this task you can use :ref:`moduleIn3Utils:Initialize Project`
* copy ``com6RockAvalanche/com6RockAvalancheCfg.ini`` to ``com6RockAvalanche/local_com6RockAvalancheCfg.ini`` and if desired change configuration settings
* if you are on a develop installation, make sure you have an updated compilation, see :ref:`advancedUsage:Update AvaFrame`
* run:
  ::

    python3 runCom6RockAvalanche.py

