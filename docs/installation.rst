Installation 
============

.. Note::
  There are two ways to install and use AvaFrame:

  Operational (GUI) 
    If you want the operational workflow, running *com1DFA* (dense flow avalanche) and *com2AB* (alpha beta)
    with standard settings from within *QGis*, head over to the :ref:`installation:Operational Installation`
    (:ref:`Deutsche Version<installation:Operationelle Installation (Deutsch)>`) section.
    Use this if you want to:

    - use the standard, well tested and calibrated setup for hazard mapping or similar
    - use QGis as frontend
    - have results in a short amount of time 
    - use the latest release 


  Advanced (Script) 
    If you want to contribute and develop AvaFrame, head over to :ref:`advancedUsage:Advanced (Script) Installation`.
    Use this if you want to:

    - work on the code itself
    - implement new features
    - change/improve existing code
    - have the latest development code. *Warning: might be unstable!*

..  Experiment **-Does not work at the moment; still under development-**
    If you want to build your own workflows and experiment with all modules,
    head over to the :ref:`installation:Experiment setup and run` section.
    Use this if you:

    - are familiar with programming in python and the terminal
    - want to build your own workflow
    - just want to adjust parameters in the configurations
    - want to use the latest release

Operational Installation 
------------------------

This is the quick start for the operational AvaFrame setup with QGis as
frontend. 

Requirements
^^^^^^^^^^^^

The prerequisites are:

* QGis: install from here: `QGis installation <https://qgis.org/en/site/forusers/download.html>`_ (we recommend
  using the latest version)

Setup AvaFrameConnector and run
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. Open QGis from your start menu and go to Plugins -> Manage and Install Plugins

#. Search for `AvaFrameConnector` and install it

#. Wait for it to finish and restart QGis

#. Access the QGis - Avaframe connector via Toolbox ->  AvaFrame -> Operational -> FullOperationalRun

#. Add the described data and run. Results will be loaded after a while
   (depending on the size of your DEM)

This tries to install the AvaFrame core package as well (since AvaFrameConnector version 1.2). If it fails, head 
over to the :ref:`operational:Manual Installation` instructions. 

For more information about the functions in the AvaFrameConnector, see :ref:`connector:AvaFrameConnector (GUI)`


Update Avaframe to a new release
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. Restart/Open QGis from your start menu and go to Plugins -> Manage and Install Plugins

#. Search for `AvaFrameConnector` and check whether it also needs updating

#. Access the QGis - Avaframe connector via Toolbox ->  AvaFrame 

#. Go to the *Admin* section, choose and run *Update*


--------------------------


Operationelle Installation (Deutsch)
-------------------------------

Dies ist der Schnellstart für die Verwendung von AvaFrame mit QGis als Frontend. 

Vorraussetzungen
^^^^^^^^^^^^^^^^

Die Vorraussetzungen sind:

* QGis: `QGis installation <https://qgis.org/de/site/forusers/download.html>`_ (Wir empfehlen die aktuellste Version zu 
  verwenden)

Setup AvaFrameConnector  
^^^^^^^^^^^^^^^^^^^^^^^

#. Öffnen Sie QGis über Ihr Startmenü und gehen Sie zu Erweiterungen -> Erweiterungen verwalten und installieren

#. Nach `AvaFrameConnector` suchen und installieren

#. Warten Sie, bis der Vorgang abgeschlossen ist und starten Sie QGis neu

#. Rufen Sie den QGis - AvaframeConnector über Verarbeitungswerkzeuge -> AvaFrame -> Operational -> FullOperationalRun auf

#. Fügen Sie die beschriebenen Daten hinzu und starten Sie den Connector. Die Ergebnisse werden nach einer Weile geladen 
   (abhängig von der Größe Ihres DEM)

Hiermit wird versucht das AvaFrame-Kernpaket zu installieren (seit AvaFrameConnector Version 1.2). Wenn dies fehlschlägt, gehen Sie 
zu: :ref:`operationalGerman:Manuelle Installation`.

Für eine Kurzzusammenfassung der Funktionen des AvaFrameConnector, siehe :ref:`connector:AvaFrameConnector (GUI)`.

Update Avaframe auf eine neue Version
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. Starten Sie QGis neu/öffnen Sie es über Ihr Startmenü und gehen Sie zu Plugins -> Plugins verwalten und installieren

#. Suchen Sie nach AvaFrameConnector und prüfen Sie, ob es aktualisiert werden muss

#. Rufen Sie den QGis - AvaframeConnector über Verarbeitungswerkzeuge -> AvaFrame auf

#. Dann -> Admin -> Update aufrufen 
