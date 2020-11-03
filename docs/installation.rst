Installation
================
This is a quick guide on how to install AvaFrame and the required dependencies on your machine.
AvaFrame is developed on **Linux machines** (Ubuntu/Manjaro/Arch) with recent
Python versions 3.8/3.9. Instructions for Windows/MacOS will follow at a later
stage. These instructions assume you are familiar with working in a terminal. 

Install `git <https://github.com/git-guides/install-git>`_ and python, we 
suggest to work with miniconda/anaconda. For installation see `miniconda
<https://docs.conda.io/en/latest/miniconda.html>`_ or
`anaconda <https://docs.anaconda.com/anaconda/install/linux/>`_.

Create a new `conda environment
<https://conda.io/projects/conda/en/latest/user-guide/concepts/environments.html>`_
for AvaFrame and activate it::

  conda create --name avaframe_env
  conda activate avaframe_env

Clone the AvaFrame repository and change into it::

  git clone https://github.com/avaframe/AvaFrame.git
  cd AvaFrame

AvaFrame needs following requirements as specified in the requirements.txt:

.. include:: ../requirements.txt
   :literal: 

Install these with (conda-forge is needed for cmocean)::
 
  conda install -c conda-forge --file requirements.txt

or alternatively::

  pip install -r requirements.txt

Finally install avaframe either with conda (might require a ``conda install conda-build``)::

  conda develop .

or by using pip::

  pip install -e .

This installs avaframe in editable mode, so every time you import avaframe the
current (local) version will be used.

Test it by starting ``python`` and do an ``import avaframe``. If no error comes
up, you are good to go. 
