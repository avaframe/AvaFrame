Data Visualisation
==========================


Main functions for creating visualisations of AvaFrame simulation results can be found in
:py:mod:`out3Plot` and :py:mod:`out1Peak`.
:py:mod:`com1DFA` also offers the possibility to export data on particle properties for visualisation using the open-source  application `**ParaView**<https://www.paraview.org/>`_.
In order to start analysing particle data of :py:mod:`com1DFA`, follow these steps:

in your local copy of ``com1DFA/com1DFACfg.ini``:

* in section [GENERAL] add `particles` to the `resType`

* in section [GENERAL] set which time steps you want to save at `tSteps`

* in section [VISUALISATION] set `writePartToCSV = True` and add the particles' properties
  you are interested in to the `particleProperties`

in ParaView

* open the particles file located in `data/avaDir/Outputs/com1DFA/particlesCSV` - there will be one
  available for each simulation performed containing the info on all saving time steps

* apply the filter **TableToPoints** and set X, Y, Z for the X, Y, Z columns

* switching to the **RenderView** you can now look at all the different particle properties you have saved
  for all exported time steps

* there is also the option to save an animation over all time steps
