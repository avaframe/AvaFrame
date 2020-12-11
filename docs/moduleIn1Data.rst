in1Data: computeFromDistribution
==========================

This function allows to retrieve a sample of values that are distributed following a uniform, or Pert distribution.
The sample ranges from the provided minimum to maximum value for which the PDF (probability density function) equals zero.
If the those should be excluded as their probability to occur is zero, set flagMinMax to False.


Inputs
-------

* configuration for the distribution


Outputs
--------

* sample derived from the chosen distribution
* various plots of the chosen distribution and the sample characteristics


To run
-------

* copy ``computeFromDistributionCfg.py`` to ``local_computeFromDistributionCfg.py`` (if not, the standard settings are used)
* adjust path to the desired ``NameOfAvalanche/`` folder in `AvaFrame/avaframe/avaframeCfg.ini``
* in ``AvaFrame/avaframe/`` run::

      python3 runComputeDist.py

.. _Theory:

Theory
-----------
