Experiment setup and run
------------------------

.. Warning::
   This section is currently under development and incomplete! It does not provide the 
   option to adjust configurations. Please use the develop installation for now. 

This is a quick guide on how to install AvaFrame
and the required dependencies on your machine. AvaFrame is developed on **Linux
machines** (Ubuntu/Manjaro/Arch) with recent Python versions 3.8/3.9/3.10.
These instructions assume you are familiar with working in a terminal. This
guide is currently described for *Linux* only, but expert users should be able
to adapt it to *Windows*.

Requirements
^^^^^^^^^^^^

Install python, we suggest to work with miniconda/anaconda. For installation see
`miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ or
`anaconda <https://docs.anaconda.com/anaconda/install/linux/>`_.

Setup AvaFrame
^^^^^^^^^^^^^^

Create a new `conda environment
<https://conda.io/projects/conda/en/latest/user-guide/concepts/environments.html>`_
for AvaFrame, activate it and install pip in this environment::

  conda create --name avaframe_env
  conda activate avaframe_env
  conda install pip

Install avaframe and its requirements::

  pip install avaframe

This installs the latest stable release of avaframe.

To install a prerelease, e.g. a release candiate, use::

  pip install avaframe --pre

Test it by starting ``python`` and do an ``import avaframe``. If no error comes
up, you are good to go.


