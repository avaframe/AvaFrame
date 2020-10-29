Installation
================
This is a quick guide on how to install AvaFrame and the required dependencies on your machine.

First make sure `git <https://github.com/git-guides/install-git>`_ is installed on your machine.
Then clone the Avaframe repository::

  git clone https://github.com/avaframe/AvaFrame.git

Now that you have the Avaframe code, you need to check that you have installed all the required dependencies.


Avaframe relies on the following python libraries:

* numpy
* matplotlib
* pyshp
* scipy
* cmocean
* seaborn
* pytest

If you are working on Linux, we suggest to work with anaconda (which you can install following `these instructions <https://docs.anaconda.com/anaconda/install/linux/>`_)
and work in a specific avaframe environment. Find out more about conda environments `here <https://conda.io/projects/conda/en/latest/user-guide/concepts/environments.html>`.

First, create a new environment for avaframe::

    conda create --name avaframe_env

Don't forget to activate it before going on::

    conda activate avaframe_env

Move to ``Avaframe/`` repository and install the different modules by using the provided installation setup file::

    conda install --name avaframe_env --file avaframe_env_spec.txt

If you are not using conda, you can perform the installation with the equivalent pip commands.

Finally install avaframe::

  conda develop .

Or if you are using pip ::

  pip install -e .

This installs avaframe in editable mode, so every time you import avaframe the
current (local) version will be used.


Build the documentation
------------------------

If you want to work on the documentation you need to install *Sphinx*. If you have followed the conda installation using
``avaframe_env_spec.txt``, you can omit the following steps. If not, you can install Sphinx, the *ReadTheDocs*
theme, and the *sphinxcontrib-bibtex*, which we use to include references, by running::

  pip install sphinx
  pip install sphinx-rtd-theme
  pip install sphinxcontrib-bibtex

In order to build the documentation you need to install make ::

  sudo apt install make

Then go to the ``docs\`` directory and run::

  make html

Html files of the documentation can be found in the _build directory.
