Advanced Installation (Windows)
------------------------------

This is a quick guide on how to install AvaFrame
and the required dependencies on your machine. AvaFrame is developed on **Linux
machines** (Ubuntu/Manjaro/Arch) with recent Python versions > 3.8. Caution: currently unavailable for Version 3.13
These instructions assume you are familiar with working in a terminal. This
guide is described for **Windows**. For *Linux*, see :ref:`developinstall:Advanced Installation (Linux)`.

Requirements
^^^^^^^^^^^^

Install `git <https://github.com/git-guides/install-git>`_ and python, we
suggest to work with miniconda/anaconda. For installation see `miniconda
<https://docs.conda.io/en/latest/miniconda.html>`_ or
`anaconda <https://docs.anaconda.com/anaconda/install/linux/>`_.

Install `Microsoft C++ compiler <https://wiki.python.org/moin/WindowsCompilers>`_.
Follow the installation steps for the version corresponding to the installed python version.

Setup AvaFrame
^^^^^^^^^^^^^^

Open conda shell and create a new `conda environment
<https://conda.io/projects/conda/en/latest/user-guide/concepts/environments.html>`_
for AvaFrame, activate it and install pip, numpy and cython in this environment::

  conda create --name avaframe_env
  conda activate avaframe_env
  conda install pip numpy cython

Clone the AvaFrame repository (in a directory of your choice: [YOURDIR]) and change into it::

  cd [YOURDIR]
  git clone https://github.com/avaframe/AvaFrame.git
  cd AvaFrame

Compile the cython com1DFA part::

   python setup.py build_ext --inplace

.. Warning::
   You will have to do this compilation every time something changes in the cython code. We also suggest
   to do this everytime updates from the repositories are pulled.

   **Before** compilation in Windows, make sure to delete ``AvaFrame/build`` directory, in addition to any .pyd, .c, and
   .pycache files in ``AvaFrame/avaframe/com1DFA``

Install avaframe and its requirements by **either** doing::

  pip install -e .

or if this fails (see `github issue 986 <https://github.com/avaframe/AvaFrame/issues/986>`_), do::

  python setup.py develop

This installs avaframe in editable mode, so every time you import avaframe the
current (local) version will be used.

Test it by starting ``python`` and do an ``import avaframe``. If no error comes
up, you are good to go.

If wanted: Integrate conda environment into editor of your choice (i.e. PyCharm, Spyder or any other)

Head over to :ref:`gettingstarted:First run` for the next steps.

Using QGIS from Conda
^^^^^^^^^^^^^^
It is possible to have the script installation and the plugin installed on your machine at the same time.
However, be aware: Depending on your mode of QGis installation (direct installer, conda, OSGeo4W...) your script AvaFrame installation might be overruled.  

The steps below install QGis and the Connector plugin in a separate *conda* environment. If the steps above have been followed, this method will result in two script installations where one will be used by QGis and the other can be worked on separately, preventing any overlap.

If you have another QGis installation somewhere, make sure that the AvaFrame plugin is not installed.
To create and activate a new environment, which will contain your installation of QGis open your conda terminal and run::

  conda create -n qgis_latest
  conda activate qgis_latest

To **install the latest QGIS version**, use::

  conda install qgis --channel conda-forge

.. Note::
   Conda will always try to install the latest version. If you want to use another version, you need to specify it, e.g.:

  ``conda install qgis=3.34.11 --channel conda-forge``

Now install the avaframe requirements::

  conda install pip numpy cython

Now find the directory where your environment is located. The path should be something like:
``C:\Users\USER\miniconda3\envs\qgis_latest``. Change into it::

  cd \miniconda3\envs\qgis_latest

Clone the AvaFrame repository and change into it::

  git clone https://github.com/avaframe/AvaFrame.git
  cd AvaFrame

.. Note::
  This will pull the latest developer version of AvaFrame, if you want the current release version for the plugin,
  you need to specify the version, e.g.:

  ``https://github.com/avaframe/AvaFrame/tree/1.8.3``

Compile the cython com1DFA part and create an editable avaframe version in your environment::

  python setup.py build_ext --inplace
  pip install -e .

Run QGIS::

  qgis

Now you can **install the AvaFrameConnector plugin** via QGIS as per usual (:ref:`installation:Operational Installation`).

AvaFrame should now be installed, and both installations should work and be editable separately. To display the installed version and installation location, use
``pip show avaframe`` in the respective environment, or use the *getVersion* function in the plugin. Note that the script
versions are separate, and changes in one installation will not affect the other.

**Updating QGIS**

To update QGIS to the most recent version, you need to run the following command with the respective environment active::

  conda update qgis -c conda-forge

