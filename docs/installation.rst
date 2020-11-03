Installation
================
This is a quick guide on how to install AvaFrame and the required dependencies on your machine.
AvaFrame is `tested  <https://codecov.io/gh/avaframe/AvaFrame>`_ with Python version 3.8 on Linux.

First make sure that `git <https://github.com/git-guides/install-git>`_ is installed.
Then clone the AvaFrame repository::

  git clone https://github.com/avaframe/AvaFrame.git

Now that you have the AvaFrame code, you need to check that you have installed all the required dependencies.


AvaFrame relies on the following python libraries:

* numpy
* matplotlib
* pyshp
* scipy
* cmocean
* seaborn
* pytest

If you are working on Linux, we suggest to work with anaconda (which you can install following `these instructions <https://docs.anaconda.com/anaconda/install/linux/>`_)
and create a specific environment to work with AvaFrame. Find out more about conda environments `here <https://conda.io/projects/conda/en/latest/user-guide/concepts/environments.html>`_.

First, create a new environment for AvaFrame::

    conda create --name avaframe_env

Don't forget to activate it before going on::

    conda activate avaframe_env

Move to ``AvaFrame/`` repository and install the different modules by using the provided installation setup file::

    conda install --name avaframe_env --file avaframe_env_spec.txt

If you are not using conda, you can perform the installation with the equivalent pip commands.

Finally install avaframe::

  conda develop .

Or if you are using pip ::

  pip install -e .

This installs avaframe in editable mode, so every time you import avaframe the
current (local) version will be used.
