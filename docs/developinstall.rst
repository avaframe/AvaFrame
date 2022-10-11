Advanced (Script) Installation 
------------------------------

This is a quick guide on how to install AvaFrame
and the required dependencies on your machine. AvaFrame is developed on **Linux
machines** (Ubuntu/Manjaro/Arch) with recent Python versions > 3.8.
These instructions assume you are familiar with working in a terminal. This
guide is currently described for *Linux* only, but expert users should be able
to adapt it to *Windows*.

Requirements
^^^^^^^^^^^^

Install `git <https://github.com/git-guides/install-git>`_ and python, we
suggest to work with miniconda/anaconda. For installation see `miniconda
<https://docs.conda.io/en/latest/miniconda.html>`_ or
`anaconda <https://docs.anaconda.com/anaconda/install/linux/>`_.
Some operating systems might require the python headers (e.g python-dev on ubuntu) or other supporting
libraries/packages (e.g. Visual Studio on Windows). 


Setup AvaFrame
^^^^^^^^^^^^^^

Create a new `conda environment
<https://conda.io/projects/conda/en/latest/user-guide/concepts/environments.html>`_
for AvaFrame, activate it and install pip, numpy and cython in this environment::

  conda create --name avaframe_env
  conda activate avaframe_env
  conda install pip numpy cython

Clone the AvaFrame repository (in a directory of your choice: [YOURDIR]) and change into it::

  cd [YOURDIR]
  git clone https://github.com/avaframe/AvaFrame.git
  cd AvaFrame

Compile the cython com1DFA part. You might also have to install a c-compiler (gcc or similar) through your systems
package manager::

  python setup.py build_ext --inplace

.. Warning::
   You will have to do this compilation every time something changes in the cython code. We also suggest
   to do this everytime updates from the repositories are pulled.

Install avaframe and its requirements::

  pip install -e .

This installs avaframe in editable mode, so every time you import avaframe the
current (local) version will be used.

Test it by starting ``python`` and do an ``import avaframe``. If no error comes
up, you are good to go.

Head over to :ref:`gettingstarted:First run` for the next steps.

Update AvaFrame
^^^^^^^^^^^^^^^

To update go to your avaframe repository [YOURDIR]/Avaframe,  pull the latest changes via::

  git pull

and repeat the compilation step from above::
  
  python setup.py build_ext --inplace

If there are updates on the requirements inside ``setup.py``, it might be necessary to run::

  pip install -e . 

again to get the additional requirements installed. 
