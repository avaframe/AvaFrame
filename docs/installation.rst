How to install
================

Get AvaFrame from github repository
------------------------------------

First make sure `git <https://github.com/git-guides/install-git>`_ is installed on your machine ::

    sudo apt-get update
    sudo apt-get install git-all

Clone repository::

  git clone https://github.com/avaframe/AvaFrame.git


Install dependencies
---------------------

Avaframe relies on the following python library:

* numpy
* matplotlib
* pyshp
* scipy
* cmocean
* seaborn
* pytest

We suggest to work with conda (which you can install following `these instructions<https://docs.anaconda.com/anaconda/install/linux/>`_) 
and create a specific avaframe `environment <https://conda.io/projects/conda/en/latest/user-guide/concepts/environments.html>`_ ::

    conda create --name avaframe_env

Don't forget to activate it before going on::

    conda activate avaframe_env

Move to ``Avaframe/`` repository and install the different modules::

    conda install --name avaframe_env --file avaframe_env_spec.txt

Or the equivalent pip commands.

Finally install avaframe::

  conda develop .

Or if you are using pip ::

  pip install -e .

This installs avaframe in editable mode, so everytime you import avaframe the
current (local) version will be used.
