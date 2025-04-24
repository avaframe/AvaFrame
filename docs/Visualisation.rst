Data Visualisation
==========================


Main functions for creating visualisations of AvaFrame simulation results can be found in
:py:mod:`out3Plot` and :py:mod:`out1Peak`.
:py:mod:`com1DFA` also offers the possibility to export data on particle properties for visualisation
using the open-source  application `ParaView <https://www.paraview.org/>`_.
In order to start analysing particle data of :py:mod:`com1DFA`, follow these steps:

in your local copy of ``com1DFA/com1DFACfg.ini``:

* *section [GENERAL]*: add `particles` to the `resType`

* *section [GENERAL]*: set which time steps you want to save at `tSteps`

* *section [VISUALISATION]*: set `writePartToCSV = True` and add the particles' properties
  you are interested in to the `particleProperties`

in ParaView

* open the particles file located in `data/avaDir/Outputs/com1DFA/particlesCSV` - there will be one
  available for each simulation performed containing the info on all saving time steps

* apply the filter **TableToPoints** and set X, Y, Z for the X, Y, Z columns

* switch to the **RenderView**:  you can look at all the different particle properties you have saved
  for all exported time steps

* it is also possible to save an animation over all time steps


Gallery
-------

This page contains example plots and pointers to the relevant plotting functions/modules. This is more to show which
type of plots exists, not to show all available plots!

**Currently work in progress**

.. list-table:: Example plots (click to view bigger version)
   :widths: 30 30 30
   :header-rows: 0

   * - .. image:: _static/gallery/avaAlr_ROC.png
                :target: _static/gallery/avaAlr_ROC.png
     - .. image:: _static/gallery/DEMplot3D.png
                :target: _static/gallery/DEMplot3D.png
     - .. image:: _static/gallery/ContourComparisonAIMEC.png
                :target: _static/gallery/ContourComparisonAIMEC.png
   * -
     - DEM plot; runscripts/runPlotTopo.py
     - Contour comparison; ana3Aimec/
   * - .. image:: _static/gallery/RunoutComparison.png
                :target: _static/gallery/RunoutComparison.png
     - .. image:: _static/gallery/VariableComparison.png
                :target: _static/gallery/VariableComparison.png
     - .. image:: _static/gallery/ThalwegTime.png
                :target: _static/gallery/ThalwegTime.png
   * - Runout line comparison
     - Variable comparison
     - Thalweg time
   * - .. image:: _static/gallery/ValongThalweg.png
                :target: _static/gallery/ValongThalweg.png
     - .. image:: _static/gallery/pftAlonProfileDEM.png
                :target: _static/gallery/pftAlonProfileDEM.png
     -
   * - Velocity along thalweg; outAIMEC -> plotVelThAlongThalweg
     - Profil along given line; out3Plot -> outProfilePlots
     -

