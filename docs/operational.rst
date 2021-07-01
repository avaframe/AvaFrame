Operational setup and run
=========================

This is the quick start for the operational AvaFrame setup with QGis as
frontend.

.. Note::  If you want to develop or run AvaFrame via python see the
           :ref:`getStarted:Quick start` page.

Requirements
------------

The prerequisites are:

* QGis: install from here: `QGis installation <https://qgis.org/en/site/forusers/download.html>`_ (we recommend
  using the latest version)
* Git: `see here <https://git-scm.com/downloads>`_

Setup AvaFrame
--------------

#. Open the command prompt -> Windows start menu -> cmd
#. Clone the necessary repositories with

    .. code-block:: bash

      git clone https://github.com/avaframe/AvaFrame.git
      git clone https://github.com/avaframe/QGisAF.git

#. Close the command prompt and open the `OSGeo4W Shell` from your start menu

#. Change to your HOME directory

    .. code-block:: bash

      cd %userprofile%

#. Install cython

    .. code-block:: bash

      python -m pip install cython

#. Change to the AvaFrame directory

    .. code-block:: bash

      cd Avaframe

#. Install the requirements

    .. code-block:: bash

      python -m pip install -r requirements.txt

#. Compile the com1DFA executable

    .. code-block::

      python avaframe\com1DFA\setup.py build_ext --inplace

#. Finally install AvaFrame

    .. code-block::

      python -m pip install -e .

Setup QGis
----------

#. Open QGis from your start menu

#. Add the QGis Connector directory (previous step 3) to your processing toolbox (TODO)

#. Reload the processing toolbox (or restart QGis)

#. Access the QGis - Avaframe connector via Toolbox -> scripts -> AvaFrame -> AvaFrame  QGis
