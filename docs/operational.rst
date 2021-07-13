Operational setup and run
-------------------------

This is the quick start for the operational AvaFrame setup with QGis as
frontend. Currently described for **Windows** only. 

Requirements
^^^^^^^^^^^^

The prerequisites are:

* QGis: install from here: `QGis installation <https://qgis.org/en/site/forusers/download.html>`_ (we recommend
  using the latest version)
* Git: `see here <https://git-scm.com/downloads>`_

Setup AvaFrame
^^^^^^^^^^^^^^

#. Open the command prompt -> Windows start menu -> cmd
#. Clone the necessary repositories with

    .. code-block:: bash

      git clone https://github.com/avaframe/AvaFrame.git
      git clone https://github.com/avaframe/QGisAF.git

#. Close the command prompt and open the **OSGeo4W Shell** from your start menu

#. Change to your HOME directory

    .. code-block:: bash

      cd %userprofile%

#. Activate the python 3 environment:

    .. code-block:: bash

      py3_env

#. Change to the AvaFrame directory

    .. code-block:: bash

      cd AvaFrame

#. Install the requirements (the `--user` avoids admin rights for installation)

    .. code-block:: bash

      python3 -m pip install --user -r requirements.txt

#. Compile the com1DFA cython executable. **If this fails with an error, see the
   additional info at** :ref:`operational:Troubleshooting cython compilation`.

    .. code-block::

      python3 avaframe\com1DFA\setup.py build_ext --inplace

#. Finally install AvaFrame

    .. code-block::

      python3 -m pip install --user -e .

#. Run the model::

     cd avaframe
     python3 runOperational.py


.. Note::
   If you encounter an error similar to this::

      ValueError: numpy.ndarray size changed, may indicate binary
      incompatibility. Expected 88 from C header, got 80 from PyObject

   run the following in OSGeo4W shell::

     python3 -m pip install --user --upgrade numpy pandas


Troubleshooting cython compilation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In case you run into problems with the cython compilation (due to a missing c++
compiler on Windows, etc), use following steps:

#. Find your python version in the **OSGeo4W Shell** with:

    .. code-block::

      python3 --version

#. Go to the latest release on the `github release
   page <https://github.com/avaframe/AvaFrame/releases>`_. Choose the
   *DFAfunctionCython.XXX* file from the assets according to your python version, i.e. for
   python 3.7 choose *cp37*, for python 3.9 *cp39* and so on.

#. Download this file and put it in the ``avaframe/com1DFA`` directory.





Setup QGis and run
^^^^^^^^^^^^^^^^^^

#. Open QGis from your start menu

#. Add the QGis Connector directory (previous step 2) to your processing toolbox 

#. Reload the processing toolbox (or restart QGis)

#. Access the QGis - Avaframe connector via Toolbox -> scripts -> AvaFrame -> AvaFrame  QGis

#. Add the described data and run. Results will be loaded after a while
   (depending on the size of your DEM).
