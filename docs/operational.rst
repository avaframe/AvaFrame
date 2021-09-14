Standard setup and run
----------------------

This is the quick start for the operational AvaFrame setup with QGis as
frontend. Currently described for **Windows** only (Linux user can use any terminal to go through the same steps below)

Requirements
^^^^^^^^^^^^

The prerequisites are:

* QGis: install from here: `QGis installation <https://qgis.org/en/site/forusers/download.html>`_ (we recommend
  using the latest version)

Setup AvaFrame
^^^^^^^^^^^^^^

#. Open the **OSGeo4W Shell** from your start menu (included in any Windows QGis installation)

#. Activate the python 3 environment:

    .. code-block:: bash

      py3_env

#. Install AvaFrame

    .. code-block::

      python3 -m pip install --user avaframe

.. Note::
   If you encounter an error similar to this::

      ValueError: numpy.ndarray size changed, may indicate binary
      incompatibility. Expected 88 from C header, got 80 from PyObject

   run the following in OSGeo4W shell::

     python3 -m pip install --user --upgrade numpy pandas



Setup QGis and run
^^^^^^^^^^^^^^^^^^

#. Open QGis from your start menu and go to Plugins -> Manage and Install Plugins

#. Search for `AvaFrameConnector` and install it (you might have to enable experimental plugins)

#. Access the QGis - Avaframe connector via Toolbox ->  AvaFrame -> AvaFrameConnector

#. Add the described data and run. Results will be loaded after a while
   (depending on the size of your DEM).
